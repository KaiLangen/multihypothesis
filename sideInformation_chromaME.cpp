#include <iostream>
#include <climits>
#include <cstring>

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

/**
 * Three Step Search Algorithm
 * description: full-pixel motion search. Looks for the block in ref
 *     that is closest to the block in trg at the center coordinates.
 * param:
 *     trg       - target frame
 *     ref       - reference frame
 *     mv        - motion vector reference object
 *     step      - starting step size for TSS
 *     center    - coordinates of the upper-left pixel in the target block
 *     width     - width of the frame(s)
 *     height    - height of the frame(s)
 *     blockSize - size of the block (height and width).
 *
 *
 */
void SideInformation::TSS(imgpel* trgU, imgpel* trgV,
                          imgpel* refU, imgpel* refV,
                          mvinfo& mv, int step, int center)
{
  // search start location
  unsigned int costs[9] = {UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,
                 UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX};
  int locations[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  int loc, og, cx, cy, x, y;
  og = center;
  // calculate the first center
  // avoid recalculating the center within the loop
  costs[4] = calcSAD(&trgU[og], &refU[center], _width, _blockSize);
  costs[4] += calcSAD(&trgV[og], &refV[center], _width, _blockSize);
  locations[4] = center;

  while(step >= 1)
  {
    // center coordinates in the image = (cy, cx)
    cy = center / _width;
    cx = center % _width;
    // coordinates in the cost matrix = (i,j)
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        // the 9 pts formed by stepping away from the center = (y, x)
        //
        // (cy-step, cx-step), (cy-step,      cx), (cy-step, cx+step)
        // (cy     , cx-step), (cy     ,      cx), (cy     , cx+step)
        // (cy+step, cx-step), (cy+step,      cx), (cy+step, cx+step)
        y = cy + (i-1) * step;
        x = cx + (j-1) * step;

        // check if the pt coordinates fall outside of the image
        if(x < 0 || x >= _width - _blockSize ||
           y < 0 || y >= _height - _blockSize ||
           (i == 1 && j == 1))
        {
          continue;
        }
        costs[i*3 + j] = calcSAD(&trgU[og],
                                 &refU[y*_width + x],
                                 _width,
                                 _blockSize);
        costs[i*3 + j] += calcSAD(&trgV[og],
                                 &refV[y*_width + x],
                                 _width,
                                 _blockSize);
        locations[i*3 + j] = y*_width + x;
      }
    }
    // re-center the search window on the local minimum
    loc = find_min(costs);
    center = locations[loc];
    step /= 2;

    // set the center and location
    costs[4] = costs[loc];
    locations[4] = center;
  }

  x = og % _width;
  y = og / _width;
  cx = center % _width;
  cy = center / _width;
  // set old coordinates in MV
  mv.iCx = x;
  mv.iCy = y;
  // MV is new location - original location
  mv.iMvx = cx - x;
  mv.iMvy = cy - y;
}

void
SideInformation::ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
                    mvinfo& mv, int p, int center, int pad)
{
  // search start location
  int cx, cy, x, y, loc;
  unsigned int cost;
  unsigned int mincost = UINT_MAX;
  int paddedWidth = (_width+2*pad);
  cy = center / paddedWidth;
  cx = center % paddedWidth;
  for(int i = -p; i < p; ++i)
  {
    for(int j = -p; j < p; ++j)
    {
      y = cy + i;
      x = cx + j;

      cost = calcSAD(&trgU[center],
                     &refU[y*paddedWidth + x],
                     paddedWidth,
                     _blockSize);
      cost += calcSAD(&trgV[center],
                      &refV[y*paddedWidth + x],
                      paddedWidth,
                      _blockSize);
      if (cost < mincost) {
        loc = y * paddedWidth + x;
        mincost = cost;
      }
    }
  }

  // set the center and location
  x = loc % paddedWidth;
  y = loc / paddedWidth;
  mv.iCx = cx;
  mv.iCy = cy;
  // MV is new location - original location
  mv.iMvx = x - cx;
  mv.iMvy = y - cy;
}

void
SideInformation::ME(imgpel* refFrameU, imgpel* currFrameU,
                    imgpel* refFrameV, imgpel* currFrameV,
                    mvinfo* candidate)
{
  int idx;
  int cnt = 0;
  int pad = 40;
  int paddedWidth = (2*pad + _width);
  for (int y = pad; y <= pad + _height - _blockSize; y += _blockSize) {
    for (int x = pad; x <= pad + _width - _blockSize; x += _blockSize) {
      idx = y * paddedWidth + x;
      ES(currFrameU, currFrameV, refFrameU, refFrameV,
         candidate[cnt], _p, idx, 40);
      cnt++;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::chroma_MEMC(RefBuffer* refFrames, imgpel* sideInfo)
{
  if (_p == 0) {
    memcpy(sideInfo, refFrames->_prevKeyFrame[0], 3*(_frameSize>>1));
    return;
  }
  int ww = _width>>1;
  int hh = _height>>1;
  int chsize = _frameSize>>2;
  int padSize = 40;
  imgpel* currUChroma = new imgpel[_frameSize];
  imgpel* currVChroma = new imgpel[_frameSize];
  imgpel* mc1 = new imgpel[_frameSize];
  imgpel* mc2 = new imgpel[_frameSize];
  auto currFrame = refFrames->_currFrame;
  auto prevKey = refFrames->_prevKeyFrame;
  auto nextKey = refFrames->_nextKeyFrame;
  imgpel* currChroma = currFrame[0] + _frameSize;

  bilinear(currChroma, currUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(currChroma+chsize, currVChroma, ww, hh, ww, hh, 0, 0);
  pad(currUChroma, currFrame[1], padSize);
  pad(currVChroma, currFrame[2], padSize);
  pad(currFrame[0], currFrame[3], padSize);

  vector<mvinfo*> mvs;
  vector<imgpel*> refs;
  // prev Key
  mvs.push_back(new mvinfo[_nmv]);
  refs.push_back(prevKey[3]);
  ME(prevKey[1], currFrame[1], prevKey[2], currFrame[2], mvs.back());
  for (int iter = 0; iter < _ss; iter++) {
    spatialSmooth(prevKey[1], prevKey[2], currFrame[1], currFrame[2],
                  mvs.back(), _blockSize, padSize); 
  }

  // next Key
  mvs.push_back(new mvinfo[_nmv]);
  refs.push_back(nextKey[3]);
  ME(nextKey[1], currFrame[1], nextKey[2], currFrame[2], mvs.back());
  for (int iter = 0; iter < _ss; iter++) {
    spatialSmooth(nextKey[1], nextKey[2], currFrame[1], currFrame[2],
                  mvs.back(), _blockSize, padSize);
  }

  // reconstructed WZ frames
  for(auto it = refFrames->begin(); it != refFrames->end(); it++)
  {
    mvs.push_back(new mvinfo[_nmv]);
    refs.push_back((*it)[3]);
    ME((*it)[1], currFrame[1], (*it)[2], currFrame[2], mvs.back());
    for (int iter = 0; iter < _ss; iter++) {
      spatialSmooth((*it)[1], (*it)[2], currFrame[1], currFrame[2],
                    mvs.back(), _blockSize, padSize); 
    }
  }
  MC(sideInfo, mvs, refs, padSize);
  MC(mc1, mvs[0], refs[0], padSize);
  MC(mc2, mvs[1], refs[1], padSize);

  _model->correlationNoiseModeling(mc1, mc2);
  delete [] currUChroma;
  delete [] currVChroma;
  for (auto m : mvs)
    delete [] m;
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

  bilinear(rU, rUBuffer, _width+2*iPadSize, _height+2*iPadSize,
                         _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(rV, rVBuffer, _width+2*iPadSize, _height+2*iPadSize,
                         _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(cU, cUBuffer, _width+2*iPadSize, _height+2*iPadSize,
                         _width+2*iPadSize, _height+2*iPadSize, 0, 0);
  bilinear(cV, cVBuffer, _width+2*iPadSize, _height+2*iPadSize,
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
SideInformation::MC(imgpel* imgDst, vector<mvinfo*> mvs,
                    vector<imgpel*> refs, int padSize)
{
  int cX, cY, mvX, mvY;
  double fWeightSum, fDist;
  double* fTmp = new double[_blockSize*_blockSize];
  mvinfo* candidate;
  int paddedWidth = (2*padSize + _width);
  imgpel* ref;
  for (int i = 0; i < _nmv; i++) {
    // init vals to zero
    fWeightSum = 0.0;
    for (int j = 0; j < _blockSize*_blockSize; j++)
        fTmp[j] = 0.0;
    
    for (size_t fIdx = 0; fIdx < mvs.size(); fIdx++) {
      fDist = mvs[fIdx]->fDist;
      fWeightSum += (1/(fDist+(float)0.001));
    }

    for (size_t fIdx = 0; fIdx < mvs.size(); fIdx++) {
      candidate = mvs[fIdx];
      cX   = candidate[i].iCx - padSize;
      cY   = candidate[i].iCy - padSize;
      mvX  = candidate[i].iCx + candidate[i].iMvx;
      mvY  = candidate[i].iCy + candidate[i].iMvy;
      fDist = mvs[fIdx]->fDist;
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
      imgDst[cX+k+(cY+j)*_width] = (imgpel)(fTmp[k+j*_blockSize] /
                                            fWeightSum);
      }
    }
  }
  delete [] fTmp;
}
#endif 

void
SideInformation::MC(imgpel* imgDst, mvinfo* candidate,
                    imgpel* ref, int padSize)
{
  int cX, cY, mvX, mvY;
  int paddedWidth = (2*padSize + _width);
  for (int i = 0; i < _nmv; i++) {
    cX   = candidate[i].iCx - padSize;
    cY   = candidate[i].iCy - padSize;
    mvX  = candidate[i].iCx + candidate[i].iMvx;
    mvY  = candidate[i].iCy + candidate[i].iMvy;
    for (int j = 0; j < _blockSize; j++) {
      for (int k = 0; k < _blockSize; k++) {
        imgDst[cX+k+(cY+j)*_width] = ref[k+mvX+(j+mvY)*paddedWidth];
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void SideInformation::initPrevNextBuffers(RefBuffer* refFrames) {
  int ww = _width>>1;
  int hh = _height>>1;
  int chsize = _frameSize>>2;
  vector<imgpel*> prev = refFrames->_prevKeyFrame;
  vector<imgpel*> next = refFrames->_nextKeyFrame;
  imgpel* prevU = prev[0] + _frameSize;
  imgpel* prevV = prevU + chsize;
  imgpel* nextU = next[0] + _frameSize;
  imgpel* nextV = nextU + chsize;
  imgpel* prevUChroma = new imgpel[_frameSize];
  imgpel* prevVChroma = new imgpel[_frameSize];
  imgpel* nextUChroma = new imgpel[_frameSize];
  imgpel* nextVChroma = new imgpel[_frameSize];

  bilinear(prevU, prevUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(prevV, prevVChroma, ww, hh, ww, hh, 0, 0);
  bilinear(nextU, nextUChroma, ww, hh, ww, hh, 0, 0);
  bilinear(nextV, nextVChroma, ww, hh, ww, hh, 0, 0);

  pad(prevUChroma, prev[1], 40);
  pad(prevVChroma, prev[2], 40);
  pad(prev[0], prev[3], 40);
  pad(nextUChroma, next[1], 40);
  pad(nextVChroma, next[2], 40);
  pad(next[0], next[3], 40);

  delete [] prevUChroma;
  delete [] nextUChroma;
  delete [] prevVChroma;
  delete [] nextVChroma;
}
