#include <string.h>
#include <climits>
#include <cstdlib>
#include <opencv2/opencv.hpp>

#include "codec.h"
#include "sideInformation.h"

using namespace cv;
class Codec;

SideInformation::SideInformation(Codec* codec, CorrModel* model)
{
  _codec     = codec; 
  _model     = model;
  _width     = _codec->getFrameWidth();
  _height    = _codec->getFrameHeight();
  _frameSize = _width * _height;
#ifdef OBMC
  _blockSize = 8;
#else
  _blockSize = 16;
#endif
  _nmv       = _width * _height / (_blockSize * _blockSize);
  _mvs       = new mvinfo[_nmv];
}

void
SideInformation::createSideInfo(imgpel* prevChroma, imgpel* currChroma,
                                imgpel* imgPrevKey, imgpel* imgCurrFrame)
{
//  memcpy(imgCurrFrame, imgPrevKey, width * height);
  /* read all the key frames */
  imgpel* refUChroma = new imgpel[_frameSize];
  imgpel* refVChroma = new imgpel[_frameSize];
  imgpel* currUChroma = new imgpel[_frameSize];
  imgpel* currVChroma = new imgpel[_frameSize];

  /* upsample the Chroma into new buffer */
  Mat mchromaU_in(_height>>1, _width>>1, CV_8UC1, currChroma);
  Mat mchromaU_out(_height, _width, CV_8UC1, currUChroma);
  resize(mchromaU_in, mchromaU_out, mchromaU_out.size(), 0, 0, INTER_CUBIC); 
  Mat mchromaV_in(_height>>1, _width>>1, CV_8UC1, currChroma + ((_frameSize)>>2));
  Mat mchromaV_out(_height, _width, CV_8UC1, currVChroma);
  resize(mchromaV_in, mchromaV_out, mchromaV_out.size(), 0, 0, INTER_CUBIC); 

  Mat mrefU_in(_height>>1, _width>>1, CV_8UC1, prevChroma);
  Mat mrefU_out(_height, _width, CV_8UC1, refUChroma);
  resize(mrefU_in, mrefU_out, mrefU_out.size(), 0, 0, INTER_CUBIC); 
  Mat mrefV_in(_height>>1, _width>>1, CV_8UC1, prevChroma + ((_frameSize)>>2));
  Mat mrefV_out(_height, _width, CV_8UC1, refVChroma);
  resize(mrefV_in, mrefV_out, mrefV_out.size(), 0, 0, INTER_CUBIC); 

  ME(refUChroma, currUChroma, refVChroma, currVChroma, _mvs);

  MC(imgPrevKey, imgCurrFrame);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcSAD(imgpel* blk1, imgpel* blk2, int width, int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+x+y*width);
      imgpel pel2=*(blk2+x+y*width);
      sad+=abs(pel1-pel2);
    }
  return sad;
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
                    mvinfo* mvs)
{
  int idx;
  unsigned int sad;
  mvinfo mv;
  int cnt = 0;
  int p = 7;
  /* for every block, search each reference frame and find the best matching block. */
  /* TODO: Would be more computationally efficient to have refs be the outter loop */
  for (int y = 0; y < _height - _blockSize + 1; y += _blockSize) {
    for (int x = 0; x < _width - _blockSize + 1; x += _blockSize) {
      sad = UINT_MAX;
      idx = y * _width + x;
      ES(currFrameU, currFrameV, refFrameU, refFrameV, mv, p, idx);
      if (mv.SAD < sad) {
        mvs[cnt] = mv;
        mvs[cnt].frameNo = 0;
        sad = mv.SAD;
      }
      cnt++;
    }
  }
}

#if RESIDUAL_CODING
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

#endif

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
SideInformation::MC(imgpel* imgPrev, imgpel* imgDst)
{
  int cX, cY, mvX[3], mvY[3];
  int cand[3];
  for (int i = 0; i < _nmv; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    cX   = _mvs[i].iCx;
    cY   = _mvs[i].iCy;

    // cand0 is the current block
    mvX[0]  = cX + _mvs[i].iMvx;
    mvY[0]  = cY + _mvs[i].iMvy;
    // cand1 is above or below
    if (cY == 0) {
      mvX[1]  = cX + _mvs[i + _width / _blockSize].iMvx;
      mvY[1]  = cY + _mvs[i + _width / _blockSize].iMvy;
    } else {
      mvX[1]  = cX + _mvs[i - _width / _blockSize].iMvx;
      mvY[1]  = cY + _mvs[i - _width / _blockSize].iMvy;
    }
    // cand2 is left or right
    if (cX % _width == 0) {
      mvX[2]  = cX + _mvs[i + 1].iMvx;
      mvY[2]  = cY + _mvs[i + 1].iMvy;
    } else {
      mvX[2]  = cX + _mvs[i - 1].iMvx;
      mvY[2]  = cY + _mvs[i - 1].iMvy;
    }

    // average all 3 candidates
    for (int j = 0; j < _blockSize; j++)
      for (int k = 0; k < _blockSize; k++) {
        for(int c = 0; c < 3; c++)
          cand[c] = _H[c][j][k] * imgPrev[mvX[c] + k + (mvY[c] + j) * _width];

        imgDst[cX + k + (cY + j) * _width] =
          (cand[0] + cand[1] + cand[2] + 4) / 8;
      }
  }
}
#else 

void
SideInformation::MC(imgpel* imgPrev, imgpel* imgDst)
{
  int cX, cY, mvX, mvY;
  for (int i = 0; i < _nmv; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // get the frame that the motion vector references, then increment it
    cX   = _mvs[i].iCx;
    cY   = _mvs[i].iCy;
    mvX  = cX + _mvs[i].iMvx;
    mvY  = cY + _mvs[i].iMvy;
    
    for (int j = 0; j < _blockSize; j++) {
      memcpy(imgDst + cX + (cY + j) * _width,
             imgPrev + mvX + (mvY + j) * _width, 
             _blockSize);
    }
  }
}
#endif 
