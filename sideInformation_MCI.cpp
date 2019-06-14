
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <vector>

#include "calculations.h"
#include "codec.h"
#include "corrModel.h"
#include "decoder.h"
#include "sideInformation.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
SideInformation::SideInformation(Codec* codec, CorrModel* model,
                                 std::map<std::string, std::string> configMap)
{
  _codec = codec;
  _model = model;
  _width = _codec->getFrameWidth();
  _height = _codec->getFrameHeight();
  _frameSize = _width * _height;
#if OBMC
  _blockSize = 8;
  cout << "Overlapped Block Motion Compensation" << endl;
#else
  _blockSize = atoi(configMap["BlockSize"].c_str());
#endif
  _p = atoi(configMap["SearchWindowSize"].c_str());
  _ss = atoi(configMap["SpatialSmoothing"].c_str());
  _nmv = _width * _height / (_blockSize * _blockSize);

  // for SI refinement
  _varList0 = new mvinfo[_width * _height / 64];
  _varList1 = new mvinfo[_width * _height / 64];
  _refinedMask = new int[_width * _height / 16];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::lowpassFilter(imgpel* src, imgpel* dst, int boxSize)
{
  const int left = boxSize/2;
  const int right = boxSize - (left + 1);

  for (int y = 0; y < _height; y++)
    for (int x = 0; x < _width; x++) {
      int sum   = 0;
      int total = 0;
  
      // Sum up pixels within the box
      for (int j = y-left; j <= y+right; j++)
        for (int i = x-left; i <= x+right; i++)
          if (i >= 0 && i < _width &&
              j >= 0 && j < _height) {
            sum += src[i+j*_width];
            total++;
          }
 
      dst[x+y*_width] = (imgpel)(sum/total);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
*get approximately true motion field
*/
void SideInformation::forwardME(imgpel* prev, imgpel* curr,
                                mvinfo* candidate, const int iRange)
{
  int px,py;
  float fDist;
  float fMinDist = -1.0;
  float fRatio;
  int iMinx=0,iMiny=0;
  int iIndex;
  int iMVx,iMVy;

  mvinfo *mv= new mvinfo[_width*_height/(iRange*iRange)];

  //half pixel buffer
  imgpel* buffer= new imgpel[iRange*iRange*4];


  for(int y=0;y<_height;y=y+iRange)
    for(int x=0;x<_width;x=x+iRange)
    {
      iIndex    = x/iRange+(y/iRange)*(_width/(iRange));
      fMinDist = -1;
      //first iteration: using 16x16 block size

      for(int pos=0;pos<(32+1)*(32+1);pos++)
      {
        iMVx=static_cast<Decoder*>(_codec)->getSpiralHpelSearchX()[pos];
        iMVy=static_cast<Decoder*>(_codec)->getSpiralHpelSearchY()[pos];

        px=iMVx+x;
        py=iMVy+y;

        fRatio=float(1.0+0.05*sqrtf((float)iMVx*iMVx+(float)iMVy*iMVy));
        if(px>=0 && px<_width && py>=0 && py<_height)
        {
          fDist=(float)calcSAD(prev,curr,px,py,x,y,iRange,_width,_height);
          fDist*=fRatio;

          if(fMinDist<0 || fDist<fMinDist)
          {
            fMinDist=fDist;
            iMinx=iMVx;
            iMiny=iMVy;
          }
        }
      }

      //half pixel ME
      //bilinear intepolation
      bilinear(prev,buffer,iRange,iRange,_width,_height,x+iMinx,y+iMiny);
      fMinDist=-1;

      for(int j=0;j<2;j++)
        for(int i=0;i<2;i++)
        {
          fDist=0;
          fDist=(float)calcSAD((buffer+i+j*(2*iRange)),
                               (curr+x+y*_width),iRange*2,
                               _width,2,1,iRange);
          if(fMinDist<0 || fDist<fMinDist)
          {
            fMinDist = fDist;
            mv[iIndex].iMvx=(2*iMinx+i)/2;
            mv[iIndex].iMvy=(2*iMiny+j)/2;
            mv[iIndex].iCx=x+iMinx/2;
            mv[iIndex].iCy=y+iMiny/2;
          }
        }

    }

  //select suitable motion vector for each block
  for(int j=0;j<_height/(iRange);j++)
    for(int i=0;i<_width/(iRange);i++)
    {
      iIndex=i+j*(_width/(iRange));
      candidate[iIndex].iCx=i*iRange;
      candidate[iIndex].iCy=j*iRange;
      candidate[iIndex].iMvx=mv[iIndex].iMvx;
      candidate[iIndex].iMvy=mv[iIndex].iMvy;

      int max_area=0;
      //select closest motion vector
      for(int k=0;k<(_width*_height/(iRange*iRange));k++)
      {
        int x=candidate[iIndex].iCx;
        int y=candidate[iIndex].iCy;
        int w=(x>mv[k].iCx)?(mv[k].iCx-x+iRange):(x-mv[k].iCx+iRange);
        int h=(y>mv[k].iCy)?(mv[k].iCy-y+iRange):(y-mv[k].iCy+iRange);

        if(w>=0 && h>=0 && w*h>=0.5*iRange*iRange){
          if(max_area==0 || w*h>max_area)
          {
            candidate[iIndex].iMvx=mv[k].iMvx;
            candidate[iIndex].iMvy=mv[k].iMvy;
          }
        }
      }
    }

  delete [] mv;
  delete [] buffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* create side information of the current frame
*/
void SideInformation::sideInfoMCI(imgpel* imgPrevKey,
                                  imgpel* imgNextKey,
                                  imgpel* imgCurrFrame)
{
  imgpel* mc1 = new imgpel[_width*_height];
  imgpel* mc2 = new imgpel[_width*_height];

  imgpel* mc1_f = new imgpel[_width*_height];
  imgpel* mc1_b = new imgpel[_width*_height];
  imgpel* mc2_f = new imgpel[_width*_height];
  imgpel* mc2_b = new imgpel[_width*_height];

  createSideInfoProcess(imgPrevKey, imgNextKey, mc1_f, mc1_b, 0);
  createSideInfoProcess(imgNextKey, imgPrevKey, mc2_b, mc2_f, 1);

  for (int iy = 0; iy < _height; iy++)
    for (int ix = 0; ix < _width; ix++) {
      int i = ix+iy*(_width);

      imgCurrFrame[i] = (mc1_f[i] + mc1_b[i] + mc2_f[i] + mc2_b[i] + 2)/4;

      mc1[i] = (mc1_f[i] + mc2_f[i] + 1)/2;
      mc2[i] = (mc1_b[i] + mc2_b[i] + 1)/2;
    }

  _model->correlationNoiseModeling(mc1, mc2);

  delete [] mc1;
  delete [] mc2;
  delete [] mc1_f;
  delete [] mc1_b;
  delete [] mc2_f;
  delete [] mc2_b;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* main process of creating side information
*/
void SideInformation::createSideInfoProcess(imgpel* imgPrevKey,
                                            imgpel* imgNextKey,
                                            imgpel* imgMCForward,
                                            imgpel* imgMCBackward,
                                            int iMode)
{
  int iRange = 16;

  imgpel *imgPrevLowPass  = new imgpel[_width*_height];
  imgpel *imgNextLowPass  = new imgpel[_width*_height];

  imgpel *imgPrevKeyPadded  = new imgpel[(_width+80)*(_height+80)];
  imgpel *imgNextKeyPadded  = new imgpel[(_width+80)*(_height+80)];

  mvinfo *varCandidate        = new mvinfo[_width*_height/(iRange*iRange)];
  mvinfo *varCandidate_iter2  = (iMode == 0) ? _varList0 : _varList1;

  //apply lowpass filter to both key frames
  lowpassFilter(imgPrevKey , imgPrevLowPass, 3);
  lowpassFilter(imgNextKey , imgNextLowPass, 3);

  forwardME(imgPrevLowPass, imgNextLowPass, varCandidate , iRange);

  pad(imgPrevKey, imgPrevKeyPadded, _width, _height, 40);
  pad(imgNextKey, imgNextKeyPadded, _width, _height, 40);

  for (int iter = 0; iter < 2; iter++)
    spatialSmooth(imgPrevKeyPadded, imgNextKeyPadded,
                  varCandidate, iRange, 40);

  //copy mv
  for (int y = 0; y < _height; y += iRange)
    for (int x = 0; x < _width; x += iRange) {
      int iIndex=(x/iRange)+(y/iRange)*(_width/iRange);

      for (int j = 0; j < iRange; j += iRange/2)
        for (int i = 0; i < iRange; i += iRange/2) {
          int index2=(x+i)/(iRange/2)+(y+j)/(iRange/2)*(_width/(iRange/2));
          varCandidate_iter2[index2].iMvx = varCandidate[iIndex].iMvx;
          varCandidate_iter2[index2].iMvy = varCandidate[iIndex].iMvy;
          varCandidate_iter2[index2].iCx  = x+i;
          varCandidate_iter2[index2].iCy  = y+j;
        }
    }

# if BIDIRECT_REFINEMENT
  bidirectME(imgPrevKeyPadded, imgNextKeyPadded,
             varCandidate_iter2, 40, iRange/2);
# endif

  for (int iter = 0; iter < 3; iter++)
    spatialSmooth(imgPrevKeyPadded, imgNextKeyPadded,
                  varCandidate_iter2, iRange/2, 40);

  MC(imgPrevKeyPadded, imgNextKeyPadded, NULL, imgMCForward, imgMCBackward,
     varCandidate_iter2, NULL, 40, iRange/2, 0);

  delete [] imgPrevLowPass;
  delete [] imgNextLowPass;
  delete [] imgPrevKeyPadded;
  delete [] imgNextKeyPadded;
  delete [] varCandidate;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
//get bidirectional motion field based on candidate mvs from previous stage (forward ME)
void SideInformation::bidirectME(imgpel* imgPrev,
                                 imgpel* imgNext,
                                 mvinfo* varCandidate,
                                 const int iPadSize,
                                 const int iRange){
  float sad,min_sad;
  float r;
  int iPosx[2]={0};
  int iPosy[2]={0};
  imgpel imgPel[2];
  int iMinx,iMiny;
  int yu,yb,xl,xr;
  int iIdx[4];
  bool bRefineFlag;

  imgpel *imgPrevBuffer,*imgNextBuffer;
  imgPrevBuffer=new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgNextBuffer=new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];

  bilinear(imgPrev, imgPrevBuffer, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(imgNext, imgNextBuffer, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize, _height+2*iPadSize, 0, 0);

  for(int iIndex=0;iIndex<_width*_height/(iRange*iRange);iIndex++)
  {
    bRefineFlag=true;
    for(int i=0;i<4;i++)
    {
      iIdx[i]=-1;
    }
    int p_x=varCandidate[iIndex].iCx/iRange;
    int p_y=varCandidate[iIndex].iCy/iRange;
    if(p_y>0)                iIdx[0] = p_x+(p_y-1)*(_width/iRange);
    if(p_x>0)                iIdx[1] = p_x-1+p_y*(_width/iRange);
    if(p_x<_width/iRange-1)  iIdx[2] = p_x+1+p_y*(_width/iRange);
    if(p_y<_height/iRange-1) iIdx[3] = p_x+(p_y+1)*(_width/iRange);

    for(int i=0;i<4;i++)
    {
      if(iIdx[i]==-1)
      {
        bRefineFlag=false;
        break;
      }
    }

    if(bRefineFlag==true)
    {
      min_sad=-1;
      yu=varCandidate[iIdx[0]].iMvy;
      yb=varCandidate[iIdx[3]].iMvy;
      xl=varCandidate[iIdx[1]].iMvx;
      xr=varCandidate[iIdx[2]].iMvx;

      iMinx=iMiny=0;
      if(yu<=yb && xl<=xr)
      {
        for(int j=yu;j<=yb;j++)
          for(int i=xl;i<=xr;i++)
          {
            iPosx[0]=2*(varCandidate[iIndex].iCx+iPadSize)+i;
            iPosy[0]=2*(varCandidate[iIndex].iCy+iPadSize)+j;
            iPosx[1]=2*(varCandidate[iIndex].iCx+iPadSize)-i;
            iPosy[1]=2*(varCandidate[iIndex].iCy+iPadSize)-j;
            sad=0;

            r=float(1.0+0.05*sqrtf(float(i*i+j*j))/2.0);
            for(int y=0;y<iRange;y++)
            {
              for(int x=0;x<iRange;x++)
              {
                imgPel[0]=imgPrevBuffer[(iPosx[0]+2*x)+
                                        (iPosy[0]+2*y)*2*(2*iPadSize+_width)];
                imgPel[1]=imgNextBuffer[(iPosx[1]+2*x)+
                                        (iPosy[1]+2*y)*2*(2*iPadSize+_width)];
                sad+=abs(imgPel[0]-imgPel[1]);
              }
            }
            sad*=r;
            if(min_sad<0 || sad<min_sad){
              iMinx=i;
              iMiny=j;
              varCandidate[iIndex].fDist=sad;
              min_sad=sad;
            }
          }
        varCandidate[iIndex].iMvx=iMinx;
        varCandidate[iIndex].iMvy=iMiny;
      }
    }
  }
  delete [] imgPrevBuffer;
  delete [] imgNextBuffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* Spatial Smoothing
*/
void SideInformation::spatialSmooth(imgpel* imgPrev,
                                    imgpel* imgNext,
                                    mvinfo* varCandidate,
                                    const int iBlockSize,
                                    const int iPadSize)
{
  int iIndex[9];
  double dWeight[9];
  int iSAD[9];
  double dMinWeight;
  int iBestMVx=0,iBestMVy=0;
  int iPx[2],iPy[2];
  int iBestIdx=0;

  imgpel *imgPrevBuffer,*imgBufferNext;
  mvinfo *varRefine = new mvinfo[_width*_height/(iBlockSize*iBlockSize)];
  imgPrevBuffer     = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgBufferNext     = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];

  bilinear(imgPrev, imgPrevBuffer, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(imgNext, imgBufferNext, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize, _height+2*iPadSize, 0, 0);

  for (int j = 0; j < _height/iBlockSize; j++)
    for (int i = 0; i < _width/iBlockSize; i++) {
      for (int k = 0; k < 9; k++)
        iIndex[k] = -1;

      if(j>0)                           
        iIndex[0]=i+(j-1)*(_width/iBlockSize);
      if(i>0)                           
        iIndex[1]=i-1+j*(_width/iBlockSize);
      if(i<_width/iBlockSize-1)         
        iIndex[2]=i+1+j*(_width/iBlockSize);
      if(j<_height/iBlockSize-1) 
        iIndex[3]=i+(j+1)*(_width/iBlockSize);
      if(i>0 && j>0)
        iIndex[5]=(i-1)+(j-1)*(_width/iBlockSize);
      if(i>0 && j<_height/iBlockSize-1)
        iIndex[6]=(i-1)+(j+1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j>0) 
        iIndex[7]=(i+1)+(j-1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j<_height/iBlockSize-1)
        iIndex[8]=(i+1)+(j+1)*(_width/iBlockSize);

      iIndex[4]=i+j*(_width/iBlockSize);

      varRefine[iIndex[4]].iMvx = varCandidate[iIndex[4]].iMvx;
      varRefine[iIndex[4]].iMvy = varCandidate[iIndex[4]].iMvy;
      varRefine[iIndex[4]].iCx  = varCandidate[iIndex[4]].iCx;
      varRefine[iIndex[4]].iCy  = varCandidate[iIndex[4]].iCy;

      for (int k = 0; k < 9; k++) {
        if (iIndex[k] != -1) {
          iPx[0]=2*(varCandidate[iIndex[4]].iCx+iPadSize)+
                    varCandidate[iIndex[k]].iMvx;
          iPy[0]=2*(varCandidate[iIndex[4]].iCy+iPadSize)+
                    varCandidate[iIndex[k]].iMvy;
          iPx[1]=2*(varCandidate[iIndex[4]].iCx+iPadSize)-
                    varCandidate[iIndex[k]].iMvx;
          iPy[1]=2*(varCandidate[iIndex[4]].iCy+iPadSize)-
                    varCandidate[iIndex[k]].iMvy;

          iSAD[k] = 0;
          iSAD[k] = calcSAD(imgPrevBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+_width),
                            imgBufferNext+iPx[1]+iPy[1]*2*(2*iPadSize+_width),
                            2*(_width+2*iPadSize), 2*(_width+2*iPadSize),
                            2, 2, iBlockSize);
        }
      }

      dMinWeight = -1;

      for (int il = 0; il < 9; il++) {
        dWeight[il] = 0;

        if (iIndex[il] != -1) {
          for (int im = 0; im < 9; im++) {
            if (il != im && iIndex[im] != -1)
              dWeight[il] += (double)(iSAD[il]/std::max<double>(iSAD[im],0.0001))
                             *(abs(varCandidate[iIndex[il]].iMvx-
                                   varCandidate[iIndex[im]].iMvx)
                              +abs(varCandidate[iIndex[il]].iMvy-
                                   varCandidate[iIndex[im]].iMvy));
          }
        }
        else
          dWeight[il]=-1;

        if ((dMinWeight<0 || dWeight[il]<=dMinWeight) && dWeight[il]>=0) {
          iBestMVx = varCandidate[iIndex[il]].iMvx;
          iBestMVy = varCandidate[iIndex[il]].iMvy;
          dMinWeight = dWeight[il];
          iBestIdx = il;
        }
      }

      varRefine[iIndex[4]].iMvx = iBestMVx;
      varRefine[iIndex[4]].iMvy = iBestMVy;
      varRefine[iIndex[4]].fDist = (float)iSAD[iBestIdx];
    }

  for (int idx = 0; idx < _width*_height/(iBlockSize*iBlockSize); idx++) {
    varCandidate[idx].iMvx = varRefine[idx].iMvx;
    varCandidate[idx].iMvy = varRefine[idx].iMvy;
    varCandidate[idx].iCx  = varRefine[idx].iCx;
    varCandidate[idx].iCy  = varRefine[idx].iCy;
    varCandidate[idx].fDist= varRefine[idx].fDist;
  }

  delete [] imgPrevBuffer;
  delete [] imgBufferNext;
  delete [] varRefine;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* this function support two modes
* mode 0 bilateral motion compensation
*        require only one motion vector set
* Param :candidate
* mode 1 bidirectional motion compensation
       require two motion vector sets (forward & backward)
* Param :candidate, candiate2
*/
void SideInformation::MC(imgpel* imgPrev, imgpel* imgNext, imgpel* imgDst,
                         imgpel* imgMCf,imgpel* imgMCb, mvinfo* varCandidate,
                         mvinfo* varCandidate2, const int iPadSize,
                         const int iRange, const int iMode){
  int px[2],py[2];
  imgpel pel[2];
  imgpel *imgPrevBuffer,*imgNextBuffer;
  imgPrevBuffer=new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgNextBuffer=new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];

  bilinear(imgPrev,imgPrevBuffer, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize,_height+2*iPadSize, 0, 0);
  bilinear(imgNext,imgNextBuffer, _width+2*iPadSize, _height+2*iPadSize,
           _width+2*iPadSize, _height+2*iPadSize, 0, 0);

  for(int idx=0;idx<_width*_height/(iRange*iRange);idx++)
  {

    for(int j=0;j<iRange;j++)
      for(int i=0;i<iRange;i++)
      {
        if(iMode==0)
        {
          px[0]=2*(varCandidate[idx].iCx+i+iPadSize)+varCandidate[idx].iMvx;
          py[0]=2*(varCandidate[idx].iCy+j+iPadSize)+varCandidate[idx].iMvy;
          px[1]=2*(varCandidate[idx].iCx+i+iPadSize)-varCandidate[idx].iMvx;
          py[1]=2*(varCandidate[idx].iCy+j+iPadSize)-varCandidate[idx].iMvy;
        }
        else
        {
          px[0]=2*(varCandidate[idx].iCx+i+iPadSize)+varCandidate[idx].iMvx;
          py[0]=2*(varCandidate[idx].iCy+j+iPadSize)+varCandidate[idx].iMvy;
          px[1]=2*(varCandidate2[idx].iCx+i+iPadSize)+varCandidate2[idx].iMvx;
          py[1]=2*(varCandidate2[idx].iCy+j+iPadSize)+varCandidate2[idx].iMvy;
        }
        pel[0]=imgPrevBuffer[px[0]+py[0]*2*(2*iPadSize+_width)];
        pel[1]=imgNextBuffer[px[1]+py[1]*2*(2*iPadSize+_width)];

        int pos=(varCandidate[idx].iCx+i)+
                (varCandidate[idx].iCy+j)*_width;
        if(imgDst!=NULL)
        {
          imgDst[pos] = (pel[0]+pel[1]+1)/2;
        }
        imgMCf[pos] = pel[0];
        imgMCb[pos] = pel[1];
      }
  }
  delete [] imgPrevBuffer;
  delete [] imgNextBuffer;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getResidualFrame(imgpel* bRef, imgpel* fRef,
                                       imgpel* curr, int* residue, int* rcList)
{
  int iBlock = 8;
  int blockCount = 0;
  imgpel* refFrame;

  for (int j = 0; j < _height; j += iBlock)
    for (int i = 0; i < _width; i += iBlock) {
      refFrame = (rcList[blockCount] == 0) ? bRef : fRef;

      for (int y = 0; y < iBlock; y++)
        for (int x = 0; x < iBlock; x++) {
          int idx = (i+x) + (j+y)*_width;

          residue[idx] = curr[idx] - refFrame[idx];
        }

      blockCount++;
    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getRecFrame(imgpel *bRef, imgpel *fRef,
                                  int *iResidue, imgpel *recFrame, int* rcList)
{
  int iBlock = 8;
  int iIndex=0;
  int iPos;
  imgpel* refFrame;
  for(int j=0;j<_height;j+=iBlock)
    for(int i=0;i<_width;i+=iBlock)
    {
      refFrame = (rcList[iBlock] == 0) ? bRef : fRef;

      for(int y=0;y<iBlock;y++)
        for(int x=0;x<iBlock;x++)
        {
          iPos = (i+x) + (j+y)*_width;
          recFrame[iPos]=Clip(0,255,iResidue[iPos]+refFrame[iPos]);
        }
      iIndex++;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getSkippedRecFrame(imgpel* imgPrevKey,
                                         imgpel * imgWZFrame, int* skipMask)
{
  for(int ij=0;ij<_height;ij+=4)
    for(int ii=0;ii<_width;ii+=4)
    {
      int idx= ii/4 + ij/4*(_width/4);
      if(skipMask[idx]==1)//skip
      {
        for(int iy=0;iy<4;iy++)
          for(int ix=0;ix<4;ix++)
          {
            imgWZFrame[(ii+ix)+(ij+iy)*_width]
              =imgPrevKey[(ii+ix)+(ij+iy)*_width];
          }
      }
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getRefinedSideInfo(imgpel *imgPrevKey,
                                         imgpel *imgNextKey,
                                         imgpel *imgCurrFrame,
                                         imgpel* imgTmpRec,
                                         imgpel *imgRefined,
                                         int iMode)
{
  int iBlock       = 8;
  int iPadSize     = 40;

  imgpel *imgForward,*imgBackward;
  imgpel *imgPrevBuffer,*imgNextBuffer,*imgPrevPadded,*imgNextPadded;
  imgForward    = new imgpel[_width*_height];
  imgBackward   = new imgpel[_width*_height];
  imgPrevPadded = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)];
  imgNextPadded = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)];

  mvinfo * mvList0 = new mvinfo[_width*_height/16];
  mvinfo * mvList1 = new mvinfo[_width*_height/16];

  pad(imgPrevKey, imgPrevPadded, _width, _height, iPadSize);
  pad(imgNextKey, imgNextPadded, _width, _height, iPadSize);

  imgPrevBuffer = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgNextBuffer = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];

  bilinear(imgPrevPadded, imgPrevBuffer, _width+2*iPadSize,
           _height+2*iPadSize, _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(imgNextPadded, imgNextBuffer, _width+2*iPadSize,
           _height+2*iPadSize, _width+2*iPadSize, _height+2*iPadSize, 0, 0);

  memset(imgForward,0x0,_width*_height);
  memset(imgBackward,0x0,_width*_height);

  for(int y=0;y<_height;y+=iBlock)
    for(int x=0;x<_width;x+=iBlock)
    {
      int iIndex=(x/iBlock)+(y/iBlock)*(_width/iBlock);
      for(int j=0;j<iBlock;j+=iBlock/2)
        for(int i=0;i<iBlock;i+=iBlock/2)
        {
          int index2=(x+i)/(iBlock/2)+(y+j)/(iBlock/2)*(_width/(iBlock/2));
          mvList0[index2].iMvx = _varList0[iIndex].iMvx;
          mvList0[index2].iMvy = _varList0[iIndex].iMvy;
          mvList0[index2].iCx  = x+i;
          mvList0[index2].iCy  = y+j;

          mvList1[index2].iMvx = _varList1[iIndex].iMvx;
          mvList1[index2].iMvy = _varList1[iIndex].iMvy;
          mvList1[index2].iCx  = x+i;
          mvList1[index2].iCy  = y+j;
        }
    }

  memset(_refinedMask,0x00,4*(_width*_height/(16)));

  getRefinedSideInfoProcess(imgPrevBuffer, imgTmpRec, imgCurrFrame,
                            imgForward, mvList0, iMode);
  getRefinedSideInfoProcess(imgNextBuffer, imgTmpRec, imgCurrFrame,
                            imgBackward, mvList1, iMode);

  for(int iy=0;iy<_height;iy++)
    for(int ix=0;ix<_width;ix++)
    {
      int iIdx = ix+iy*(_width);
      imgRefined[iIdx]=(imgForward[iIdx]+imgBackward[iIdx]+1)/2;
    }
  _model->updateCNM(imgForward,imgBackward,_refinedMask);

  delete [] mvList0;
  delete [] mvList1;
  delete [] imgForward;
  delete [] imgBackward;
  delete [] imgPrevPadded;
  delete [] imgNextPadded;
  delete [] imgPrevBuffer;
  delete [] imgNextBuffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
float getDCValue(imgpel* img,int iStep,int iStripe,int iBlock)
{
  float fDCValue;
  int sum = 0;
  for(int j=0;j<iBlock;j++)
    for(int i=0;i<iBlock;i++)
    {
      sum += *(img + i*iStep + j*iStep*iStripe);
    }
  fDCValue = (float)sum/float(iBlock*iBlock);
  return fDCValue;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getRefinedSideInfoProcess(imgpel* imgPrevBuffer,
                                                imgpel* imgTmpRec,
                                                imgpel* imgSI,
                                                imgpel* imgRefined,
                                                mvinfo* varList,
                                                int iMode)
{
  int iBlock       = 4;
  int x,y,i,j;
  int iCx,iCy;
  int iPox,iPoy;
  int iPadSize     = 40;
  int iSearchRange = 8;
  int iStripe      = _width+2*iPadSize;
  float fDist      = 0.0;
  float fSIDist    = 0.0;
  float fThreshold = 3.0;
  float fRefinedTH;
  float fRefinedTH_dc = (float)2.0*_codec->getQuantStep(0, 0, false);
  float fRefinedTH_ac = 100.0;

  fRefinedTH = (!iMode)?fRefinedTH_dc:fRefinedTH_ac;
  fThreshold = (!iMode)?(float)4.0:fRefinedTH_ac;

  vector< pair<int,int> > varMVList;
  vector<float>         varWeight;
  float                 fWeightSum;
  float                 fTmp;

  //do ME again to find the better MV
  for(int iIndex=0;iIndex<(_height*_width)/(iBlock*iBlock);iIndex++)
    {
      varMVList.clear();
      varWeight.clear();
      fWeightSum = 0.0;
      x = 2*(varList[iIndex].iCx+iPadSize)+varList[iIndex].iMvx;
      y = 2*(varList[iIndex].iCy+iPadSize)+varList[iIndex].iMvy;
      iCx = varList[iIndex].iCx;
      iCy = varList[iIndex].iCy;

      float fRecDC = getDCValue(imgTmpRec+iCx+iCy*_width,1,_width,iBlock);
      float fDCValue;
      if(iMode==0)//DC
      {
        fDCValue = getDCValue(imgSI+iCx+iCy*_width,1,_width,iBlock);
        fSIDist = fabs(fDCValue - fRecDC);
      }
      else
      {
        fSIDist = (float)calcDist((imgSI+iCx+iCy*_width),
                                  (imgTmpRec+iCx+iCy*_width),
                                  _width, _width, 1, 1, iBlock);
      }
      if(fSIDist>fRefinedTH)
      {
        for(int pos=0;pos<(2*iSearchRange+1)*(2*iSearchRange+1);pos++)
        {
            i=static_cast<Decoder*>(_codec)->getSpiralSearchX()[pos];
            j=static_cast<Decoder*>(_codec)->getSpiralSearchY()[pos];

            if(iMode == 0)//DC
            {
              fDCValue = getDCValue(imgPrevBuffer+(x+i)+(y+j)*iStripe*2,
                                    2, iStripe*2, iBlock);
              fDist  =  fabs(fDCValue - fRecDC);
            }
            else
            {
              fDist = (float)calcDist((imgPrevBuffer+(x+i)+(y+j)*iStripe*2),
                                      (imgTmpRec+iCx+iCy*_width),iStripe*2,
                                      _width,2,1,iBlock);
            }

            if(pos==0 && iMode!=0) fThreshold = fDist;

            if( fDist < 0.8*fThreshold)
            {
              varMVList.push_back(pair<int,int>(i,j));
              varWeight.push_back(1/(fDist+(float)0.001));
              fWeightSum += (1/(fDist+(float)0.001));
            }
        }
      }
      //weighted mean
      for(int ij=0;ij<iBlock;ij++)
        for(int ii=0;ii<iBlock;ii++)
        {
          fTmp = 0.0;
          if(varMVList.size()!=0)
          {
            for(size_t iList = 0;iList<varMVList.size();iList++)
            {
              iPox = x + 2*ii + varMVList[iList].first ;
              iPoy = y + 2*ij + varMVList[iList].second;
              fTmp += (*(imgPrevBuffer+iPox+iPoy*iStripe*2))*varWeight[iList];
            }
            imgRefined[(iCx+ii)+(iCy+ij)*_width] = (imgpel) (fTmp/fWeightSum);
          }
          else
          {
            imgRefined[(iCx+ii)+(iCy+ij)*_width] =
                 imgSI[(iCx+ii)+(iCy+ij)*_width];
          }
        }

      if(varMVList.size()!=0)
      {
        _refinedMask[iIndex]++;
      }
    }
}
