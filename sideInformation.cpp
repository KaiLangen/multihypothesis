#include <string.h>
#include <climits>
#include <iostream>
#include <cstdlib>

#include "codec.h"
#include "decoder.h"
#include "sideInformation.h"

class Codec;


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
#else
  _blockSize = atoi(configMap["BlockSize"].c_str());
#endif
  _p = atoi(configMap["SearchWindowSize"].c_str());
  _ss = atoi(configMap["SpatialSmoothing"].c_str());
  _nmv = _width * _height / (_blockSize * _blockSize);
  _mvs = new mvinfo[_nmv];
}

void
SideInformation::createSideInfo(imgpel* prevChroma, imgpel* currChroma,
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
  bilinear(prevChroma, refUChroma, ww, hh, ww, hh);
  bilinear(prevChroma+(_frameSize>>2), refVChroma, ww, hh, ww, hh);

  bilinear(currChroma, currUChroma, ww, hh, ww, hh);
  bilinear(currChroma+(_frameSize>>2), currVChroma, ww, hh, ww, hh);
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
  int iIndex[9];
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
                        _width+2*iPadSize, _height+2*iPadSize);
  bilinear(rV,rVBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize);
  bilinear(cU,cUBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize);
  bilinear(cV,cVBuffer, _width+2*iPadSize, _height+2*iPadSize,
                        _width+2*iPadSize, _height+2*iPadSize);

  for (int j = 0; j < _height/iBlockSize; j++)
    for (int i = 0; i < _width/iBlockSize; i++) {
      for (int k = 0; k < 9; k++)
        iIndex[k] = -1;

      if(j>0)                           iIndex[0]=i+(j-1)*(_width/iBlockSize);
      if(i>0)                           iIndex[1]=i-1+j*(_width/iBlockSize);
      if(i<_width/iBlockSize-1)         iIndex[2]=i+1+j*(_width/iBlockSize);
      if(j<_height/iBlockSize-1)        iIndex[3]=i+(j+1)*(_width/iBlockSize);
      if(i>0 && j>0)                    iIndex[5]=(i-1)+(j-1)*(_width/iBlockSize);
      if(i>0 && j<_height/iBlockSize-1) iIndex[6]=(i-1)+(j+1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j>0)  iIndex[7]=(i+1)+(j-1)*(_width/iBlockSize);
      if(i<_width/iBlockSize-1 && j<_height/iBlockSize-1) iIndex[8]=(i+1)+(j+1)*(_width/iBlockSize);

      iIndex[4]=i+j*(_width/iBlockSize);

      varRefine[iIndex[4]].iMvx = varCandidate[iIndex[4]].iMvx;
      varRefine[iIndex[4]].iMvy = varCandidate[iIndex[4]].iMvy;
      varRefine[iIndex[4]].iCx  = varCandidate[iIndex[4]].iCx;
      varRefine[iIndex[4]].iCy  = varCandidate[iIndex[4]].iCy;

      for (int k = 0; k < 9; k++) {
        if (iIndex[k] != -1) {
          iPx[0]=2*(varCandidate[iIndex[4]].iCx+iPadSize)+varCandidate[iIndex[k]].iMvx;
          iPy[0]=2*(varCandidate[iIndex[4]].iCy+iPadSize)+varCandidate[iIndex[k]].iMvy;
          iPx[1]=2*(varCandidate[iIndex[4]].iCx+iPadSize)-varCandidate[iIndex[k]].iMvx;
          iPy[1]=2*(varCandidate[iIndex[4]].iCy+iPadSize)-varCandidate[iIndex[k]].iMvy;

          iSAD[k] = 0;
          iSAD[k] = calcSAD(rUBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+_width),
                            cUBuffer+iPx[1]+iPy[1]*2*(2*iPadSize+_width),
                            2*(_width+2*iPadSize), 2*(_width+2*iPadSize), 2, 2, iBlockSize);
          iSAD[k] += calcSAD(rVBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+_width),
                            cVBuffer+iPx[1]+iPy[1]*2*(2*iPadSize+_width),
                            2*(_width+2*iPadSize), 2*(_width+2*iPadSize), 2, 2, iBlockSize);
        }
      }

      dMinWeight = -1;

      for (int il = 0; il < 9; il++) {
        dWeight[il] = 0;

        if (iIndex[il] != -1) {
          for (int im = 0; im < 9; im++) {
            if (il != im && iIndex[im] != -1)
              dWeight[il] += (double)(iSAD[il]/std::max<double>(iSAD[im],0.0001))
                               *( abs(varCandidate[iIndex[il]].iMvx-varCandidate[iIndex[im]].iMvx)
                                 +abs(varCandidate[iIndex[il]].iMvy-varCandidate[iIndex[im]].iMvy));
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

  for (int iIndex = 0; iIndex < _width*_height/(iBlockSize*iBlockSize); iIndex++) {
    varCandidate[iIndex].iMvx = varRefine[iIndex].iMvx;
    varCandidate[iIndex].iMvy = varRefine[iIndex].iMvy;
    varCandidate[iIndex].iCx  = varRefine[iIndex].iCx;
    varCandidate[iIndex].iCy  = varRefine[iIndex].iCy;
    varCandidate[iIndex].fDist= varRefine[iIndex].fDist;
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
void SideInformation::getResidualFrame(imgpel* bRefFrame, imgpel* currFrame, int* residue)
{
  int iBlock = 8;
  int blockCount = 0;

  for (int j = 0; j < _height; j += iBlock)
    for (int i = 0; i < _width; i += iBlock) {
      for (int y = 0; y < iBlock; y++)
        for (int x = 0; x < iBlock; x++) {
          int idx = (i+x) + (j+y)*_width;

          residue[idx] = currFrame[idx] - bRefFrame[idx];
        }

      blockCount++;
    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getRecFrame(imgpel *imgBReference, int *iResidue, imgpel *imgRec)
{
  int iBlock = 8;
  int iIndex=0;
  int iPos;
  for(int j=0;j<_height;j+=iBlock)
    for(int i=0;i<_width;i+=iBlock) {
      for(int y=0;y<iBlock;y++)
        for(int x=0;x<iBlock;x++)
        {
          iPos = (i+x) + (j+y)*_width;
          imgRec[iPos]=Clip(0,255,iResidue[iPos]+imgBReference[iPos]);
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


void
SideInformation::MC(imgpel* imgPrev, imgpel* imgDst, int padSize)
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
              int picwidth, int picheight){
    for(int j=0;j<buffer_h;j++)
    for(int i=0;i<buffer_w;i++)
    {
      int buffer_r=2*buffer_w;
      int a,b,c,d;

      int x=i;
      int y=j;
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
