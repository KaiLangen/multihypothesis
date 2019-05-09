#include <string.h>
#include <iostream>
#include <climits>
#include <cstdlib>
#include <math.h>

#include "codec.h"
#include "decoder.h"
#include "sideInformation.h"

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

int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2, int s1,int s2,int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+s1*x+s1*y*width1);
      imgpel pel2=*(blk2+s2*x+s2*y*width2);
      sad+=abs(pel1-pel2);
    }
  return sad;
}


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

bool SideInformation::isSkip(const int* skipMask, int start, int thresh) {
  int skipCount = 0;
  int wSkipOffset = _width / 4;
  for (int y = 0; y < _blockSize/4; y++)
    for (int x = 0; x < _blockSize/4; x++)
      skipCount += skipMask[start+x+y*wSkipOffset];

  if (skipCount > thresh) return true;
  else return false;

}

void
SideInformation::ME(imgpel* refFrameU, imgpel* currFrameU,
                    imgpel* refFrameV, imgpel* currFrameV)
{
  int idx;
  int cnt = 0;
  int wSkipOffset = _width / 4;
  int* skipMask = static_cast<Decoder*>(_codec)->_skipMask;
  /* for every block, search each reference frame and find the best matching block. */
  /* TODO: Would be more computationally efficient to have refs be the outter loop */
  for (int y = 0; y <= _height - _blockSize; y += _blockSize) {
    for (int x = 0; x <= _width - _blockSize; x += _blockSize) {
      /* check skipmask subblocks */
      int skipStart = (y*4)*wSkipOffset + (x*4);
      if (isSkip(skipMask, skipStart, _blockSize/4+1)) {
        _mvs[cnt].iCx = x;
        _mvs[cnt].iCy = y;
        _mvs[cnt].iMvx = 0;
        _mvs[cnt].iMvy = 0;
      } else {
        idx = y * _width + x;
        ES(currFrameU, currFrameV, refFrameU, refFrameV, _mvs[cnt], _p, idx);
      }
      cnt++;
    }
  }
}
