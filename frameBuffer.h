#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include <vector>
#include <cstring>
#include <cstdio>

#include "defs.h"

using namespace std;
class RefBuffer
{
public:
  // Each frame is represented as 4 buffers to reduce computation
  // 1 - original YUV frame
  // 2 - Upsampled/Padded U channel
  // 3 - Upsampled/Padded V channel
  // 4 - Padded Y channel
  vector<vector<imgpel*>> _refFrames;
  vector<imgpel*> _currFrame;
  int _windowSize;
  int _nextRec;
  bool _isFull;
  int _padSize;
  int _width;
  int _height;
  int _frameSize;
  int _yuvFrameSize;
  int _paddedFrameSize;

  RefBuffer(int width, int height, int windowSize);

  ~RefBuffer();

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void updateRecWindow(FILE* fKeyReadPtr, int prevKeyNo);

  void init(FILE* fKeyReadPtr);

  vector<vector<imgpel*>>::iterator begin() { return _refFrames.begin(); }

  vector<vector<imgpel*>>::iterator end()
  {
    if (_isFull) return _refFrames.end();
    else return _refFrames.begin() + _nextRec;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
};

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int windowSize = 0);

  ~FrameBuffer();

  //  Luma
  imgpel*  getCurrFrame()        { return _rBuff->_currFrame[0]; };
  imgpel*  getorigFrame()        { return _origFrame; };
  imgpel*  getSideInfoFrame()    { return _sideInfoFrame; };
  RefBuffer* getRefBuffer()      { return _rBuff; };

private:
  int      _frameSize;
  imgpel*  _origFrame;
  imgpel*  _sideInfoFrame;
  RefBuffer* _rBuff;
};

// general functions
void pad(imgpel* src, imgpel* dst, int width, int height, int iPadSize);

void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight, int px, int py);

void decimate2(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
               int picwidth, int px, int py);
#endif // ENCODER_INC_FRAMEBUFFER_H

