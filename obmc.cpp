#include <iostream>

#include "types.h"
#include "obme.h"

/* the hardest part to figure out... */
void obmc(imgpel** refFrame, imgpel* currFrame, mvinfo* mvs,
          int width, int height, int blocksize)
{
  int idx, cX, cY, mvX, mvY, numMV;
  imgpel* ref;
  numMV = (width * height) / (blocksize * blocksize);
  for (int i = 0; i < numMV; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // get the frame that the motion vector references, then increment it
    cX  = mvs[i].iCx;
    mvX = cX + mvs[i].iMvx;
    cY  = mvs[i].iCy;
    mvY = cY + mvs[i].iMvy;
    idx = mvs[i].frame;
    ref = refFrame[idx];
    for (int y = 0; y < blocksize; y++) {
      for (int x = 0; x < blocksize; x++) {
        currFrame[(cX + x) + (cY + y) * width] += 
        ref[(mvX + x) + (mvY + y) * width];
      }
    }
  }
}
