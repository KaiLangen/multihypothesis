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
 * return: the final SAD of the matching block.
 *
 */
int SideInformation::TSS(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
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
//                    costs[i*3 + j] = USHRT_MAX;
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

    return costs[4];
}

void
SideInformation::ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
                    mvinfo& mv, int p, int center)
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
            cost = calcSAD(&trgU[center],
                           &refU[y*_width + x],
                           _width,
                           _blockSize);
            cost += calcSAD(&trgV[center],
                            &refV[y*_width + x],
                            _width,
                            _blockSize);
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
    mv.SAD = mincost;
}

void
SideInformation::ME(imgpel* refFrameU, imgpel* currFrameU,
                    imgpel* refFrameV, imgpel* currFrameV,
                    mvinfo* candidate)
{
  int idx;
  int cnt = 0;
  /* for every block, search each reference frame and find the best matching block. */
  /* TODO: Would be more computationally efficient to have refs be the outter loop */
  for (int y = 0; y <= _height - _blockSize; y += _blockSize) {
    for (int x = 0; x <= _width - _blockSize; x += _blockSize) {
      idx = y * _width + x;
      ES(currFrameU, currFrameV, refFrameU, refFrameV, candidate[cnt], _p, idx);
      cnt++;
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

  ME(prevUChroma, currUChroma, prevVChroma, currVChroma, varCandidate0);
  ME(nextUChroma, currUChroma, nextVChroma, currVChroma, varCandidate1);
  
  for (int iter = 0; iter < _ss; iter++) { 
    spatialSmooth(pUPadded, pVPadded, cUPadded, cVPadded,
                  varCandidate0, _blockSize, 40); 
    spatialSmooth(nUPadded, nVPadded, cUPadded, cVPadded,
                  varCandidate1, _blockSize, 40); 
  }

  MC(prevPadded, mcF, varCandidate0, 40);
  MC(nextPadded, mcB, varCandidate1, 40);

  for (int iy = 0; iy < _height; iy++)
    for (int ix = 0; ix < _width; ix++) {
      int i = ix+iy*(_width);
      imgCurrFrame[i] = (mcF[i] + mcB[i] + 1)/2;
    }

  _model->correlationNoiseModeling(mcF, mcB);
  
  // delete all buffers
  delete []  prevUChroma;
  delete []  prevVChroma;
  delete []  nextUChroma;
  delete []  nextVChroma;
  delete []  currUChroma;
  delete []  currVChroma;
  delete []  pUPadded;
  delete []  pVPadded;
  delete []  nUPadded;
  delete []  nVPadded;
  delete []  cUPadded;
  delete []  cVPadded;
  delete []  nextPadded;
  delete []  prevPadded;
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
