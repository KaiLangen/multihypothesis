#include <string.h>
#include <iostream>
#include <climits>

#include "calculations.h"
#include "codec.h"
#include "corrModel.h"
#include "decoder.h"
#include "sideInformation.h"


static int find_min(unsigned int costs[9])
{
    unsigned int minimum = costs[0];
    int location = 0;
    for(int c = 1; c < 9; ++c)
    {
        if(costs[c] < minimum)
        {
            minimum = costs[c];
            location = c;
        }
    }
    return location;
}

void
SideInformation::ES(mvinfo& mv, int p, int center,
                    imgpel* ref1, imgpel* trg1,
                    imgpel* ref2=0, imgpel* trg2=0)
{
    // search start location
    int cx, cy, x, y, loc;
    unsigned int cost;
    unsigned int mincost = UINT_MAX;
    cy = center / _width;
    cx = center % _width;
    for(int i = -p; i < p; ++i)
    {
        for(int j = -p; j < p; ++j)
        {
            y = cy + i;
            x = cx + j;

            // check if the pt coordinates fall outside of the image
            if(x < 0 || x >= _width - _blockSize ||
               y < 0 || y >= _height - _blockSize ||
               (i == 1 && j == 1))
            {
                continue;
            }
            cost = calcSAD(&trg1[center],
                           &ref1[y*_width + x],
                           _width,
                           _blockSize);
            if (trg2) {
              cost += calcSAD(&trg2[center],
                              &ref2[y*_width + x],
                              _width,
                              _blockSize);
            }
            if (cost < mincost) {
              loc = y * _width + x;
              mincost = cost;
            }
        }
    }

    // set the center and location
    x = loc % _width;
    y = loc / _width;
    mv.iCx = cx;
    mv.iCy = cy;
    // MV is new location - original location
    mv.iMvx = x - cx;
    mv.iMvy = y - cy;
    mv.fDist = (float)mincost;
}

void
SideInformation::ME(imgpel* refFrameU, imgpel* currFrameU,
                    imgpel* refFrameV, imgpel* currFrameV,
                    mvinfo* candidate)
{
  int idx;
  int c = 0;
  /* for every block, search each reference
   * frame and find the best matching block. */
  for (int y = 0; y <= _height - _blockSize; y += _blockSize) {
    for (int x = 0; x <= _width - _blockSize; x += _blockSize) {
      idx = y * _width + x;
      ES(candidate[c], _p, idx, refFrameU, currFrameU, refFrameV, currFrameV);
      c++;
    }
  }
}

void
SideInformation::ME(imgpel* refFrame, imgpel* currFrame,
                    mvinfo* candidate)
{
  int idx;
  int c = 0;
  /* for every block, search each reference
   * frame and find the best matching block. */
  for (int y = 0; y <= _height - _blockSize; y += _blockSize) {
    for (int x = 0; x <= _width - _blockSize; x += _blockSize) {
      idx = y * _width + x;
      ES(candidate[c], _p, idx, refFrame, currFrame);
      c++;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::chroma_MEMC(imgpel* prevChroma, imgpel* imgPrevKey,
                                  imgpel* nextChroma, imgpel* imgNextKey,
                                  RecFrameBuffer* recFrames,
				  imgpel* currChroma, imgpel* imgCurrFrame)
{
  if (_p == 0) {
    memcpy(imgCurrFrame, imgPrevKey, _frameSize);
    return;
  }
  mvinfo *varCandidate0 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  mvinfo *varCandidate1 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  mvinfo *varCandidate2 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  mvinfo *varCandidate3 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  imgpel* prevUChroma = new imgpel[_frameSize];
  imgpel* prevVChroma = new imgpel[_frameSize];
  imgpel* nextUChroma = new imgpel[_frameSize];
  imgpel* nextVChroma = new imgpel[_frameSize];
  imgpel* currUChroma = new imgpel[_frameSize];
  imgpel* currVChroma = new imgpel[_frameSize];
  imgpel* pUPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* pVPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* nUPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* nVPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* cUPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* cVPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* nextPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* prevPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* mcF = new imgpel[_frameSize];
  imgpel* mcB = new imgpel[_frameSize];

  vector< pair<int,int> > varMVList;
  vector<float>         varWeight;

  /* upsample the Chroma into new buffer, and pad */
  int ww = _width>>1;
  int hh = _height>>1;
  int chsize = _frameSize>>2;
  bilinear(prevChroma, prevUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(prevChroma+chsize, prevVChroma, ww, hh, ww, hh, 0, 0);

  bilinear(nextChroma, nextUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(nextChroma+chsize, nextVChroma, ww, hh, ww, hh, 0, 0);

  bilinear(currChroma, currUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(currChroma+chsize, currVChroma, ww, hh, ww, hh, 0, 0);

  pad(currUChroma, cUPadded, 40);
  pad(currVChroma, cVPadded, 40);
  pad(prevUChroma, pUPadded, 40);
  pad(prevVChroma, pVPadded, 40);
  pad(nextUChroma, nUPadded, 40);
  pad(nextVChroma, nVPadded, 40);
  pad(imgPrevKey, prevPadded, 40);
  pad(imgNextKey, nextPadded, 40);

  if (isDoubleMV) {
    ME(prevUChroma, currUChroma, varCandidate0);
    ME(prevVChroma, currVChroma, varCandidate1);
    ME(nextUChroma, currUChroma, varCandidate2);
    ME(nextVChroma, currVChroma, varCandidate3);

    for (int iter = 0; iter < _ss; iter++) {
      spatialSmooth(pUPadded, pVPadded, cUPadded, cVPadded,
                    varCandidate0, _blockSize, 40);
      spatialSmooth(pUPadded, pVPadded, cUPadded, cVPadded,
                    varCandidate1, _blockSize, 40);
      spatialSmooth(nUPadded, nVPadded, cUPadded, cVPadded,
                    varCandidate2, _blockSize, 40);
      spatialSmooth(nUPadded, nVPadded, cUPadded, cVPadded,
                    varCandidate3, _blockSize, 40);
    }
    mvinfo* mvs[4] = {varCandidate0, varCandidate1, varCandidate2, varCandidate3};
    imgpel* refs[4] = {prevPadded, prevPadded, nextPadded, nextPadded};
    MC(refs, mvs, imgCurrFrame, 40, 4);
    MC(prevPadded, mcF, varCandidate0, 40);
    MC(nextPadded, mcB, varCandidate2, 40);
  } else {
    ME(prevUChroma, currUChroma, prevVChroma, currVChroma, varCandidate0);
    ME(nextUChroma, currUChroma, nextVChroma, currVChroma, varCandidate1);

    for (int iter = 0; iter < _ss; iter++) {
      spatialSmooth(pUPadded, pVPadded, cUPadded, cVPadded,
                    varCandidate0, _blockSize, 40);
      spatialSmooth(nUPadded, nVPadded, cUPadded, cVPadded,
                    varCandidate1, _blockSize, 40);
    }
    mvinfo* mvs[2] = {varCandidate0, varCandidate1};
    imgpel* refs[2] = {prevPadded, nextPadded};
    MC(refs, mvs, imgCurrFrame, 40, 2);
    MC(prevPadded, mcF, varCandidate0, 40);
    MC(nextPadded, mcB, varCandidate1, 40);
  }

  _model->correlationNoiseModeling(mcF, mcB);

  // delete all buffers
  delete [] varCandidate0;
  delete [] varCandidate1;
  delete [] varCandidate2;
  delete [] varCandidate3;
  delete [] prevUChroma;
  delete [] prevVChroma;
  delete [] nextUChroma;
  delete [] nextVChroma;
  delete [] currUChroma;
  delete [] currVChroma;
  delete [] pUPadded;
  delete [] pVPadded;
  delete [] nUPadded;
  delete [] nVPadded;
  delete [] cUPadded;
  delete [] cVPadded;
  delete [] nextPadded;
  delete [] prevPadded;
  delete [] mcF;
  delete [] mcB;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::oracle_MEMC(imgpel* imgPrevKey, imgpel* imgNextKey,
                                  imgpel* currLuma, imgpel* imgResult)
{
  if (_p == 0) {
    memcpy(imgResult, imgPrevKey, _frameSize);
    return;
  }
  mvinfo *varCandidate0 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  mvinfo *varCandidate1 = new mvinfo[_width*_height/(_blockSize*_blockSize)];
  imgpel* nextPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* prevPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* currPadded = new imgpel[(_width+80)*(_height+80)];
  imgpel* mcF = new imgpel[_frameSize];
  imgpel* mcB = new imgpel[_frameSize];

  pad(imgPrevKey, prevPadded, 40);
  pad(imgNextKey, nextPadded, 40);
  pad(currLuma, currPadded, 40);

  ME(imgPrevKey, currLuma, varCandidate0);
  ME(imgNextKey, currLuma, varCandidate1);

  for (int iter = 0; iter < _ss; iter++) {
    spatialSmooth(prevPadded, currPadded,
                  varCandidate0, _blockSize, 40);
    spatialSmooth(nextPadded, currPadded,
                  varCandidate1, _blockSize, 40);
  }

  MC(prevPadded, mcF, varCandidate0, 40);
  MC(nextPadded, mcB, varCandidate1, 40);

  mvinfo* mvs[2] = {varCandidate0, varCandidate1};
  imgpel* refs[2] = {prevPadded, nextPadded};
  MC(refs, mvs, imgResult, 40, 2);

  _model->correlationNoiseModeling(mcF, mcB);

  // delete all buffers
  delete []  nextPadded;
  delete []  prevPadded;
  delete []  currPadded;
  delete []  mcF;
  delete []  mcB;
}

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

void SideInformation::MC(imgpel* imgPrev, imgpel* imgDst,
                         mvinfo* candidate, int padSize)
{

  int cX, cY, mvX[3], mvY[3];
  int cand[3];
  for (int i = 0; i < _nmv; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    cX   = candidate[i].iCx;
    cY   = candidate[i].iCy;

    // cand0 is the current block
    mvX[0]  = cX + candidate[i].iMvx + padSize;
    mvY[0]  = cY + candidate[i].iMvy + padSize;
    // cand1 is above or below
    if (cY == 0) {
      mvX[1]  = cX + candidate[i + _width / _blockSize].iMvx + padSize;
      mvY[1]  = cY + candidate[i + _width / _blockSize].iMvy + padSize;
    } else {
      mvX[1]  = cX + candidate[i - _width / _blockSize].iMvx + padSize;
      mvY[1]  = cY + candidate[i - _width / _blockSize].iMvy + padSize;
    }
    // cand2 is left or right
    if (cX == 0) {
      mvX[2]  = cX + candidate[i + 1].iMvx + padSize;
      mvY[2]  = cY + candidate[i + 1].iMvy + padSize;
    } else {
      mvX[2]  = cX + candidate[i - 1].iMvx + padSize;
      mvY[2]  = cY + candidate[i - 1].iMvy + padSize;
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
    cX   = candidate[i].iCx;
    cY   = candidate[i].iCy;
    mvX  = cX + candidate[i].iMvx + padSize;
    mvY  = cY + candidate[i].iMvy + padSize;

    for (int j = 0; j < _blockSize; j++) {
      memcpy(imgDst + cX + (cY + j) * _width,
             imgPrev + mvX + (mvY + j) * (2*padSize+_width),
             _blockSize);
    }
  }
}
#endif

void
SideInformation::MC(imgpel* refs[], mvinfo* mvs[],
                    imgpel* imgDst, int padSize, int nRefs)
{
  int cX, cY, mvX, mvY;
  double fWeightSum, fDist;
  double* fTmp = new double[_blockSize*_blockSize];
  int paddedWidth = (2*padSize + _width);
  mvinfo* candidate;
  imgpel* ref;

  for (int i = 0; i < _nmv; i++) {
    // init vals to zero
    fWeightSum = 0.0;
    for (int j = 0; j < _blockSize*_blockSize; j++)
        fTmp[j] = 0.0;

    for (int fIdx = 0; fIdx < nRefs; fIdx++) {
      fDist = mvs[fIdx][i].fDist;
      fWeightSum += (1/(fDist+(float)0.001));
    }

    for (int fIdx = 0; fIdx < nRefs; fIdx++) {
      candidate = mvs[fIdx];
      cX   = candidate[i].iCx;
      cY   = candidate[i].iCy;
      mvX  = cX + candidate[i].iMvx + padSize;
      mvY  = cY + candidate[i].iMvy + padSize;
      fDist = candidate[i].fDist;
      ref = refs[fIdx];
      for (int j = 0; j < _blockSize; j++) {
        for (int k = 0; k < _blockSize; k++) {
          fTmp[k+j*_blockSize] += ref[k+mvX+(j+mvY)*paddedWidth] *
                                  (1/(fDist+(float)0.001));
        }
      }
    }
    for (int j = 0; j < _blockSize; j++) {
      for (int k = 0; k < _blockSize; k++) {
      imgDst[cX+k+(cY+j)*_width] = (imgpel)((fTmp[k+j*_blockSize])/
                                            fWeightSum);
      }
    }
  }
  delete [] fTmp;
}

