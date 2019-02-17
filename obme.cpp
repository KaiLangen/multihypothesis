#include <climits>
#include <cstdlib>
#include <iostream>

#include "types.h"
#include "obme.h"

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

void ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV, mvinfo& mv,
       int p, int center, int width, int height, int blocksize)
{
    // search start location
    int cx, cy, x, y, loc;
    unsigned int cost;
    unsigned int mincost = UINT_MAX;
    cy = center / width;
    cx = center % width;
    for(int i = -p; i < p; ++i)
    {
        for(int j = -p; j < p; ++j)
        {
            y = cy + i;
            x = cx + j;

            // check if the pt coordinates fall outside of the image
            if(x < 0 || x >= width - blocksize ||
               y < 0 || y >= height - blocksize ||
               (i == 1 && j == 1))
            {
                continue;
            }
            cost = calcSAD(&trgU[center],
                           &refU[y*width + x],
                           width,
                           blocksize);
            cost += calcSAD(&trgV[center],
                            &refV[y*width + x],
                            width,
                            blocksize);
            if (cost < mincost) {
              loc = y * width + x;
              mincost = cost;
            }
        }
    }

    // set the center and location
    x = loc % width;
    y = loc / width;
    mv.iCx = cx;
    mv.iCy = cy;
    // MV is new location - original location
    mv.iMvx = x - cx;
    mv.iMvy = y - cy;
    mv.SAD = mincost;
}

void obme(imgpel** refFrameU, imgpel* currFrameU,
          imgpel** refFrameV, imgpel* currFrameV, mvinfo* mvs,
          int nRefs, int width, int height, int blocksize)
{
  int idx;
  unsigned int sad;
  mvinfo mv;
  int cnt = 0;
  int p = 7;
  /* for every block, search each reference frame and find the best matching block. */
  /* TODO: Would be more computationally efficient to have refs be the outter loop */
  for (int x = 0; x < width - blocksize + 1; x += blocksize) {
    for (int y = 0; y < height - blocksize + 1; y += blocksize) {
      sad = UINT_MAX;
      for (int i = 0; i < nRefs; i++) {
        idx = y * width + x;
        ES(currFrameU, currFrameV, refFrameU[i], refFrameV[i],
           mv, p, idx, width, height, blocksize);
        if (mv.SAD < sad) {
//          std::cout << sad << " > " << mv.SAD << std::endl;
          mvs[cnt] = mv;
          mvs[cnt].frame = i;
          sad = mv.SAD;
        }
      }
      cnt++;
    }
  }
}
