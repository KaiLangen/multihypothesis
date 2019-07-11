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
                    mvinfo& mv, int p, int center, int padSize)
{
  // search start location
  int cx, cy, x, y, loc;
  unsigned int cost;
  unsigned int mincost = UINT_MAX;
  int paddedWidth = (_width+2*padSize);
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

  // padSize must be subtracted from coordinates to work with spatial smoothing
  mv.iCx = cx - padSize;
  mv.iCy = cy - padSize;
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
void SideInformation::oracle_MEMC(RefBuffer* refFrames, imgpel* sideInfo)
{
  int chsize = (_width*_height)>>2;
  int padSize = 40;
  auto currFrame = refFrames->_currFrame;

  pad(currFrame[0], currFrame[3], _width, _height, padSize);

  vector<mvinfo*> mvs;
  vector<imgpel*> refsU;
  vector<imgpel*> refsV;
  // all reference Key-frames
  for(auto it = refFrames->begin(); it != refFrames->end(); it++)
  {
    mvs.push_back(new mvinfo[_nmv]);
    refsU.push_back((*it)[1]);
    refsV.push_back((*it)[2]);
    ME((*it)[3], currFrame[3], (*it)[3], currFrame[3], mvs.back());
    for (int iter = 0; iter < _ss; iter++) {
      spatialSmooth((*it)[3], (*it)[3], currFrame[3], currFrame[3],
                    mvs.back(), _blockSize, padSize); 
    }
  }

  MC(sideInfo, mvs, refsU, padSize/2);
  MC(sideInfo+chsize, mvs, refsV, padSize/2);

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

void
SideInformation::MC(imgpel* imgDst, vector<mvinfo*> mvs,
                    vector<imgpel*> refs, int padSize)
{
  int cX, cY, mvX, mvY;
  double fWeightSum, fDist;
  int hblock = _blockSize>>1;
  int ww = _width>>1;
  double* fTmp = new double[hblock*hblock];
  mvinfo* candidate;
  int paddedWidth = (2*padSize + ww);
  imgpel* ref;
  for (int i = 0; i < _nmv; i++) {
    // init vals to zero
    fWeightSum = 0.0;
    for (int j = 0; j < hblock*hblock; j++)
        fTmp[j] = 0.0;
    
    for (size_t fIdx = 0; fIdx < mvs.size(); fIdx++) {
      fDist = mvs[fIdx][i].fDist;
      fWeightSum += (1/(fDist+(float)0.001));
    }

    for (size_t fIdx = 0; fIdx < mvs.size(); fIdx++) {
      candidate = mvs[fIdx];
      cX   = candidate[i].iCx/2;
      cY   = candidate[i].iCy/2;
      mvX  = cX + candidate[i].iMvx/2 + padSize;
      mvY  = cY + candidate[i].iMvy/2 + padSize;
      fDist = candidate[i].fDist;
      ref = refs[fIdx];
      for (int j = 0; j < hblock; j++) {
        for (int k = 0; k < hblock; k++) {
          fTmp[k+j*hblock] += ref[k+mvX+(j+mvY)*paddedWidth] *
                                  (1/(fDist+(float)0.001));
        }
      }
    }
    for (int j = 0; j < hblock; j++) {
      for (int k = 0; k < hblock; k++) {
      imgDst[cX+k+(cY+j)*ww] = (imgpel)((fTmp[k+j*hblock]) /
                                            fWeightSum);
      }
    }
  }
  delete [] fTmp;
}
