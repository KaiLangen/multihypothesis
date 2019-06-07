
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
  _mvs = new mvinfo[_nmv];

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
void SideInformation::createSideInfo(imgpel* prevChroma, imgpel* currChroma,
                                     imgpel* imgPrevKey, imgpel* imgCurrFrame)
{
  if (_p == 0) {
    memcpy(imgCurrFrame, imgPrevKey, _frameSize);
    return;
  }
  /* upsample the Chroma into new buffer */
  imgpel* refUChroma = new imgpel[_frameSize];
  imgpel* refVChroma = new imgpel[_frameSize];
  imgpel* currUChroma = new imgpel[_frameSize];
  imgpel* currVChroma = new imgpel[_frameSize];
  imgpel* rUPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* rVPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* cUPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* cVPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* prevPadded = new imgpel[(_width+80)*(_height+80)];

  int ww = _width>>1;
  int hh = _height>>1;
  bilinear(prevChroma, refUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(prevChroma+(_frameSize>>2), refVChroma, ww, hh, ww, hh, 0, 0);

  bilinear(currChroma, currUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(currChroma+(_frameSize>>2), currVChroma, ww, hh, ww, hh, 0, 0);
  pad(currUChroma, cUPadded, 40);
  pad(currVChroma, cVPadded, 40);
  pad(refUChroma, rUPadded, 40);
  pad(refVChroma, rVPadded, 40);
  pad(imgPrevKey, prevPadded, 40);

  ME(refUChroma, currUChroma, refVChroma, currVChroma);
  
  for (int iter = 0; iter < _ss; iter++) 
    spatialSmooth(rUPadded, rVPadded, cUPadded, cVPadded, _mvs, _blockSize, 40); 

  MC(prevPadded, imgCurrFrame, 40);
}
//#endif

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* Spatial Smoothing
*/
void SideInformation::spatialSmooth(imgpel* rU, imgpel* rV, imgpel* cU, imgpel*cV,
    mvinfo* varCandidate, const int iBlockSize, const int iPadSize)
{
  int ind[9];
  double dWeight[9];
  int iSAD[9];
  double dMinWeight;
  int iBestMVx=0,iBestMVy=0;
  int iPx[2],iPy[2];
  int iBestIdx=0;

  mvinfo *varRefine = new mvinfo[_width*_height/(iBlockSize*iBlockSize)];
  imgpel* rUBuffer  = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgpel* rVBuffer  = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgpel* cUBuffer  = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];
  imgpel* cVBuffer  = new imgpel[(_width+2*iPadSize)*(_height+2*iPadSize)*4];

  bilinear(rU,rUBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(rV,rVBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(cU,cUBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(cV,cVBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize, 0, 0);

  for (int j = 0; j < _height/iBlockSize; j++)
    for (int i = 0; i < _width/iBlockSize; i++) {
      for (int k = 0; k < 9; k++)
        ind[k] = -1;

      if(j>0)
        ind[0]=i+(j-1)*(_width/iBlockSize);
      if(i>0)
        ind[1]=i-1+j*(_width/iBlockSize);
      if(i<_width/iBlockSize-1)
        ind[2]=i+1+j*(_width/iBlockSize);
      if(j<_height/iBlockSize-1)
        ind[3]=i+(j+1)*(_width/iBlockSize);
      if(i>0 && j>0)
        ind[5]=(i-1)+(j-1)*(_width/iBlockSize);
      if(i>0 && j<_height/iBlockSize-1)
        ind[6]=(i-1)+(j+1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j>0)
        ind[7]=(i+1)+(j-1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j<_height/iBlockSize-1)
        ind[8]=(i+1)+(j+1)*(_width/iBlockSize);

      ind[4]=i+j*(_width/iBlockSize);

      varRefine[ind[4]].iMvx = varCandidate[ind[4]].iMvx;
      varRefine[ind[4]].iMvy = varCandidate[ind[4]].iMvy;
      varRefine[ind[4]].iCx  = varCandidate[ind[4]].iCx;
      varRefine[ind[4]].iCy  = varCandidate[ind[4]].iCy;

      for (int k = 0; k < 9; k++) {
        if (ind[k] != -1) {
          iPx[0]=2*(varCandidate[ind[4]].iCx+iPadSize)
                   +varCandidate[ind[k]].iMvx;
          iPy[0]=2*(varCandidate[ind[4]].iCy+iPadSize)
                   +varCandidate[ind[k]].iMvy;
          iPx[1]=2*(varCandidate[ind[4]].iCx+iPadSize)
                   -varCandidate[ind[k]].iMvx;
          iPy[1]=2*(varCandidate[ind[4]].iCy+iPadSize)
                   -varCandidate[ind[k]].iMvy;

          iSAD[k] = 0;
          iSAD[k] = calcSAD(rUBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+_width),
                            cUBuffer+iPx[1]+iPy[1]*2*(2*iPadSize+_width),
                            2*(_width+2*iPadSize), 2*(_width+2*iPadSize),
                            2, 2, iBlockSize);
          iSAD[k] += calcSAD(rVBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+_width),
                             cVBuffer+iPx[1]+iPy[1]*2*(2*iPadSize+_width),
                             2*(_width+2*iPadSize), 2*(_width+2*iPadSize),
                             2, 2, iBlockSize);
        }
      }

      dMinWeight = -1;

      for (int il = 0; il < 9; il++) {
        dWeight[il] = 0;

        if (ind[il] != -1) {
          for (int im = 0; im < 9; im++) {
            if (il != im && ind[im] != -1)
              dWeight[il] += (double)(iSAD[il]/
                              std::max<double>(iSAD[im],0.0001))
                              *(abs(varCandidate[ind[il]].iMvx
                                   -varCandidate[ind[im]].iMvx)
                               +abs(varCandidate[ind[il]].iMvy
                                   -varCandidate[ind[im]].iMvy));
          }
        }
        else
          dWeight[il]=-1;

        if ((dMinWeight<0 || dWeight[il]<=dMinWeight) && dWeight[il]>=0) {
          iBestMVx = varCandidate[ind[il]].iMvx;
          iBestMVy = varCandidate[ind[il]].iMvy;
          dMinWeight = dWeight[il];
          iBestIdx = il;
        }
      }

      varRefine[ind[4]].iMvx = iBestMVx;
      varRefine[ind[4]].iMvy = iBestMVy;
      varRefine[ind[4]].fDist = (float)iSAD[iBestIdx];
    }

  for (int idx = 0; idx < _width*_height/(iBlockSize*iBlockSize); idx++) {
    varCandidate[idx].iMvx = varRefine[idx].iMvx;
    varCandidate[idx].iMvy = varRefine[idx].iMvy;
    varCandidate[idx].iCx  = varRefine[idx].iCx;
    varCandidate[idx].iCy  = varRefine[idx].iCy;
    varCandidate[idx].fDist= varRefine[idx].fDist;
  }

  delete [] rUBuffer;
  delete [] rVBuffer;
  delete [] cUBuffer;
  delete [] cVBuffer;
  delete [] varRefine;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::pad(imgpel* src, imgpel* dst, const int iPadSize)
{
  int padded_width  = _width  + 2*iPadSize;
  int padded_height = _height + 2*iPadSize;

  // Loops start from iPadSize; subtract it from src index

  // Upper left
  for (int y = 0; y < iPadSize; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padded_width] = src[0];

  // Upper
  for (int x = iPadSize; x < iPadSize+_width; x++)
    for (int y = 0; y < iPadSize; y++)
      dst[x+y*padded_width] = src[x-iPadSize];

  // Upper right
  for (int y = 0; y < iPadSize; y++)
    for (int x = iPadSize+_width; x < padded_width; x++)
      dst[x+y*padded_width] = src[_width-1];

  // Left
  for (int y = iPadSize; y < iPadSize+_height; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padded_width] = src[(y-iPadSize)*_width];

  // Middle
  for (int y = iPadSize; y < iPadSize+_height; y++)
    for (int x = iPadSize; x < iPadSize+_width; x++)
      dst[x+y*padded_width] = src[x-iPadSize+(y-iPadSize)*_width];

  // Right
  for (int y = iPadSize; y < iPadSize+_height; y++)
    for (int x = iPadSize+_width; x < padded_width; x++)
      dst[x+y*padded_width] = src[(y-iPadSize+1)*_width-1];

  // Bottom left
  for (int y = iPadSize+_height; y < padded_height; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padded_width] = src[(_height-1)*_width];

  // Bottom
  for (int y = iPadSize+_height; y < padded_height; y++)
    for (int x = iPadSize; x < iPadSize+_width; x++)
      dst[x+y*padded_width] = src[(x-iPadSize)+(_height-1)*_width];

  // Bottom right
  for (int y = iPadSize+_height; y < padded_height; y++)
    for (int x = iPadSize+_width; x < padded_width; x++)
      dst[x+y*padded_width] = src[_height*_width-1];
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

#if OBMC
const int SideInformation::_H[3][8][8] =
{
  {
    {4, 5, 5, 5, 5, 5, 5, 4},
    {5, 5, 5, 5, 5, 5, 5, 5},
    {5, 5, 6, 6, 6, 6, 5, 5},
    {5, 5, 6, 6, 6, 6, 5, 5},
    {5, 5, 6, 6, 6, 6, 5, 5},
    {5, 5, 6, 6, 6, 6, 5, 5},
    {5, 5, 5, 5, 5, 5, 5, 5},
    {4, 5, 5, 5, 5, 5, 5, 4}
  },
  {
    {2, 2, 2, 2, 2, 2, 2, 2},
    {1, 1, 2, 2, 2, 2, 1, 1},
    {1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 2, 2, 2, 2, 1, 1},
    {2, 2, 2, 2, 2, 2, 2, 2}
  },
  {
    {2, 1, 1, 1, 1, 1, 1, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 2, 1, 1, 1, 1, 2, 2},
    {2, 1, 1, 1, 1, 1, 1, 2}
  }
};

void SideInformation::MC(imgpel* imgPrev, imgpel* imgDst, int padSize)
{
  
  int cX, cY, mvX[3], mvY[3];
  int cand[3];
  for (int i = 0; i < _nmv; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    cX   = _mvs[i].iCx;
    cY   = _mvs[i].iCy;

    // cand0 is the current block
    mvX[0]  = cX + _mvs[i].iMvx + padSize;
    mvY[0]  = cY + _mvs[i].iMvy + padSize;
    // cand1 is above or below
    if (cY == 0) {
      mvX[1]  = cX + _mvs[i + _width / _blockSize].iMvx + padSize;
      mvY[1]  = cY + _mvs[i + _width / _blockSize].iMvy + padSize;
    } else {
      mvX[1]  = cX + _mvs[i - _width / _blockSize].iMvx + padSize;
      mvY[1]  = cY + _mvs[i - _width / _blockSize].iMvy + padSize;
    }
    // cand2 is left or right
    if (cX == 0) {
      mvX[2]  = cX + _mvs[i + 1].iMvx + padSize;
      mvY[2]  = cY + _mvs[i + 1].iMvy + padSize;
    } else {
      mvX[2]  = cX + _mvs[i - 1].iMvx + padSize;
      mvY[2]  = cY + _mvs[i - 1].iMvy + padSize;
    }

    // average all 3 candidates
    for (int j = 0; j < _blockSize; j++)
      for (int k = 0; k < _blockSize; k++) {
        for(int c = 0; c < 3; c++)
          cand[c] = _H[c][j][k] * imgPrev[mvX[c]+k+
                                         (mvY[c]+j)*(2*padSize+_width)];

        imgDst[cX + k + (cY + j) * _width] =
          (cand[0] + cand[1] + cand[2] + 4) / 8;
      }
  }
}
#else 

void
SideInformation::MC(imgpel* imgPrev, imgpel* imgDst, int padSize)
{
  int cX, cY, mvX, mvY;
  for (int i = 0; i < _nmv; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // get the frame that the motion vector references, then increment it
    cX   = _mvs[i].iCx;
    cY   = _mvs[i].iCy;
    mvX  = cX + _mvs[i].iMvx + padSize;
    mvY  = cY + _mvs[i].iMvy + padSize;
    
    for (int j = 0; j < _blockSize; j++) {
      memcpy(imgDst + cX + (cY + j) * _width,
             imgPrev + mvX + (mvY + j) * (2*padSize+_width),
             _blockSize);
    }
  }
}
#endif 

void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight, int px, int py){
    for(int j=0;j<buffer_h;j++)
    for(int i=0;i<buffer_w;i++)
    {
      int buffer_r=2*buffer_w;
      int a,b,c,d;

      int x=px+i;
      int y=py+j;
      if(x>picwidth-1)x=picwidth-1;
      if(x<0)x=0;
      if(y>picheight-1)y=picheight-1;
      if(y<0)y=0;
      a=source[(x)+(y)*picwidth];

      if((x+1)<picwidth) b=source[(x+1)+(y)*picwidth];
      else b=a;

      if((y+1)<picheight) c=source[(x)+(y+1)*picwidth];
      else c=a;

      if((x+1)<picwidth && (y+1)<picheight) d=source[(x+1)+(y+1)*picwidth];
      else d=a;

      buffer[2*i+(2*j)*buffer_r]=a;
      buffer[(2*i+1)+(2*j)*buffer_r]=(a+b)/2;
      buffer[2*i+(2*j+1)*buffer_r]=(a+c)/2;
      buffer[(2*i+1)+(2*j+1)*buffer_r]=(a+b+c+d)/4;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getSkippedRecFrame(imgpel* imgPrevKey,
                                         imgpel * imgWZFrame, int* skipMask)
{
  int iWidth,iHeight;
  iWidth  = _codec->getFrameWidth();
  iHeight = _codec->getFrameHeight();

  for(int ij=0;ij<iHeight;ij+=4)
    for(int ii=0;ii<iWidth;ii+=4)
    {
      int idx= ii/4 + ij/4*(iWidth/4);
      if(skipMask[idx]==1)//skip
      {
        for(int iy=0;iy<4;iy++)
          for(int ix=0;ix<4;ix++)
          {
            imgWZFrame[(ii+ix)+(ij+iy)*iWidth]
              =imgPrevKey[(ii+ix)+(ij+iy)*iWidth];
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
  int iWidth,iHeight;
  iWidth           = _codec->getFrameWidth();
  iHeight          = _codec->getFrameHeight();
  int iBlock       = 8;
  int iPadSize     = 40;

  imgpel *imgForward,*imgBackward;
  imgpel *imgPrevBuffer,*imgNextBuffer,*imgPrevPadded,*imgNextPadded;
  imgForward    = new imgpel[iWidth*iHeight];
  imgBackward   = new imgpel[iWidth*iHeight];
  imgPrevPadded = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)];
  imgNextPadded = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)];

  mvinfo * mvList0 = new mvinfo[iWidth*iHeight/16];
  mvinfo * mvList1 = new mvinfo[iWidth*iHeight/16];

  pad(imgPrevKey,imgPrevPadded,iPadSize);
  pad(imgNextKey,imgNextPadded,iPadSize);

  imgPrevBuffer = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)*4];
  imgNextBuffer = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)*4];

  bilinear(imgPrevPadded, imgPrevBuffer, iWidth+2*iPadSize,
           iHeight+2*iPadSize, iWidth+2*iPadSize, iHeight+2*iPadSize, 0, 0);
  bilinear(imgNextPadded, imgNextBuffer, iWidth+2*iPadSize,
           iHeight+2*iPadSize, iWidth+2*iPadSize, iHeight+2*iPadSize, 0, 0);

  memset(imgForward,0x0,iWidth*iHeight);
  memset(imgBackward,0x0,iWidth*iHeight);

  for(int y=0;y<iHeight;y+=iBlock)
    for(int x=0;x<iWidth;x+=iBlock)
    {
      int iIndex=(x/iBlock)+(y/iBlock)*(iWidth/iBlock);
      for(int j=0;j<iBlock;j+=iBlock/2)
        for(int i=0;i<iBlock;i+=iBlock/2)
        {
          int index2=(x+i)/(iBlock/2)+(y+j)/(iBlock/2)*(iWidth/(iBlock/2));
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

  memset(_refinedMask,0x00,4*(iWidth*iHeight/(16)));

  getRefinedSideInfoProcess(imgPrevBuffer, imgTmpRec, imgCurrFrame,
                            imgForward, mvList0, iMode);
  getRefinedSideInfoProcess(imgNextBuffer, imgTmpRec, imgCurrFrame,
                            imgBackward, mvList1, iMode);

  for(int iy=0;iy<iHeight;iy++)
    for(int ix=0;ix<iWidth;ix++)
    {
      int iIdx = ix+iy*(iWidth);
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
  int iWidth,iHeight;
  iWidth           = _codec->getFrameWidth();
  iHeight          = _codec->getFrameHeight();
  int iBlock       = 4;
  int x,y,i,j;
  int iCx,iCy;
  int iPox,iPoy;
  int iPadSize     = 40;
  int iSearchRange = 8;
  int iStripe      = iWidth+2*iPadSize;
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
  for(int iIndex=0;iIndex<(iHeight*iWidth)/(iBlock*iBlock);iIndex++)
    {
      varMVList.clear();
      varWeight.clear();
      fWeightSum = 0.0;
      x = 2*(varList[iIndex].iCx+iPadSize)+varList[iIndex].iMvx;
      y = 2*(varList[iIndex].iCy+iPadSize)+varList[iIndex].iMvy;
      iCx = varList[iIndex].iCx;
      iCy = varList[iIndex].iCy;

      float fRecDC = getDCValue(imgTmpRec+iCx+iCy*iWidth,1,iWidth,iBlock);
      float fDCValue;
      if(iMode==0)//DC
      {
        fDCValue = getDCValue(imgSI+iCx+iCy*iWidth,1,iWidth,iBlock);
        fSIDist = fabs(fDCValue - fRecDC);
      }
      else
      {
        fSIDist = (float)calcDist((imgSI+iCx+iCy*iWidth),
                                  (imgTmpRec+iCx+iCy*iWidth),
                                  iWidth, iWidth, 1, 1, iBlock);
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
                                      (imgTmpRec+iCx+iCy*iWidth),iStripe*2,
                                      iWidth,2,1,iBlock);
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
            imgRefined[(iCx+ii)+(iCy+ij)*iWidth] = (imgpel) (fTmp/fWeightSum);
          }
          else
          {
            imgRefined[(iCx+ii)+(iCy+ij)*iWidth] =
                 imgSI[(iCx+ii)+(iCy+ij)*iWidth];
          }
        }

      if(varMVList.size()!=0)
      {
        _refinedMask[iIndex]++;
      }
    }
}
