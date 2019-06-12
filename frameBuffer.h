
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include <vector>

#include "defs.h"

class RecFrameBuffer
{
public:
  RecFrameBuffer(int frameSize, int windowSize)
  {
    _windowSize = windowSize;
    _currRec    = 0;
    if (windowSize > 0)
      _recFrames(windowSize, new imgpel[frameSize]);
  }

  ~RecFrameBuffer()
  {
    for (auto p : _recFrames)
      delete [] p;
  }

  // increment currRec, and wrap if equal to buffSize
  imgpel* getNextRec()
  {
    if (++_currRec >= _windowSize) _currRec = 0;
    return _recFrames[_currRec];
  }

private:
  vector<imgpel*> _recFrames;
  int _windowSize;
  int _currRef;
};

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int windowSize = 0)
  {
    (void)gop;
    _frameSize        = width * height;
    _prevKeyFrame     = new imgpel[3*(_frameSize>>1)];
    _currFrame        = new imgpel[3*(_frameSize>>1)];
    _nextKeyFrame     = new imgpel[3*(_frameSize>>1)];
    _origFrame        = new imgpel[3*(_frameSize>>1)];
    _sideInfoFrame    = new imgpel[3*(_frameSize>>1)];
    _dctFrame         = new int[3*(_frameSize>>1)];
    _quantDctFrame    = new int[3*(_frameSize>>1)];
    _decFrame         = new int[3*(_frameSize>>1)];
    _invQuantDecFrame = new int[3*(_frameSize>>1)];
    _recFrames        = new RecFrameBuffer(_frameSize, windowSize);

  };

  ~FrameBuffer()
  {
    delete [] _prevKeyFrame;
    delete [] _currFrame;
    delete [] _nextKeyFrame;
    delete [] _origFrame;
    delete [] _sideInfoFrame;
    delete [] _dctFrame;
    delete [] _quantDctFrame;
    delete [] _decFrame;
    delete [] _invQuantDecFrame;
    delete _recFrames;

  }

  //  Luma
  imgpel*  getPrevFrame()        { return _prevFrame; };
  imgpel*  getCurrFrame()        { return _currFrame; };
  imgpel*  getNextFrame()        { return _nextFrame; };
  imgpel*  getorigFrame()        { return _origFrame; };
  imgpel*  getSideInfoFrame()    { return _sideInfoFrame; };
  int*     getDctFrame()         { return _dctFrame; };
  int*     getQuantDctFrame()    { return _quantDctFrame; };
  int*     getDecFrame()         { return _decFrame; };
  int*     getInvQuantDecFrame() { return _invQuantDecFrame; };
  RecFrameBuffer* getRecFrameBuffer() { return _recFrames; };

  // Chroma
  imgpel*  getPrevChroma()
  { return _prevFrame + _frameSize; };
  imgpel*  getCurrChroma()
  { return _currFrame + _frameSize; };
  imgpel*  getNextChroma()
  { return _nextFrame + _frameSize; };
  imgpel*  getorigChroma()
  { return _origFrame + _frameSize; };
  int*     getDctChroma()
  { return _dctFrame + _frameSize; };
  int*     getQuantDctChroma()
  { return _quantDctFrame + _frameSize; };
  int*     getDecChroma()
  { return _decFrame + _frameSize; };
  int*     getInvQuantDecChroma()
  { return _invQuantDecFrame + _frameSize; };

private:
  int      _frameSize;
  imgpel*  _prevKeyFrame;
  imgpel*  _currFrame;
  imgpel*  _nextKeyFrame;
  imgpel*  _origFrame;
  imgpel*  _sideInfoFrame;

  int*     _dctFrame;
  int*     _quantDctFrame;
  int*     _decFrame;
  int*     _invQuantDecFrame;
  RecFrameBuffer* _recFrames;
};

#endif // ENCODER_INC_FRAMEBUFFER_H

