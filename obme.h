#include "types.h"

#ifndef OBME_H_
#define OBME_H_

#include "obme.h"

void lowpassFilter(imgpel* src, imgpel* dst,
                   const int boxSize, const int width, const int height);

void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight);

void obme(imgpel** refFrameU, imgpel* currFrameU,
          imgpel** refFrameV, imgpel* currFrameV, mvinfo* mvs,
          int nRefs, int width, int height, int blocksize);

void obmc(imgpel** refFrame, imgpel* currFrame, mvinfo* mvs,
          int width, int height, int blocksize);

#endif
