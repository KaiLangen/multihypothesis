
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include "defs.h"

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int gop = 0)
  {
    (void)gop;
    _frameSize = width * height;
    _prevFrame        = new imgpel[3*(_frameSize>>1)];
    _currFrame        = new imgpel[3*(_frameSize>>1)];
    _nextFrame        = new imgpel[3*(_frameSize>>1)];
    _origFrame        = new imgpel[3*(_frameSize>>1)];
    _sideInfoFrame    = new imgpel[3*(_frameSize>>1)];
    _dctFrame         = new int[3*(_frameSize>>1)];
    _quantDctFrame    = new int[3*(_frameSize>>1)];
    _decFrame         = new int[3*(_frameSize>>1)];
    _invQuantDecFrame = new int[3*(_frameSize>>1)];
  };

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
  int _frameSize;
  imgpel*  _prevFrame;
  imgpel*  _currFrame;
  imgpel*  _nextFrame;
  imgpel*  _origFrame;
  imgpel*  _sideInfoFrame;

  int*     _dctFrame;
  int*     _quantDctFrame;
  int*     _decFrame;
  int*     _invQuantDecFrame;
};

#endif // ENCODER_INC_FRAMEBUFFER_H

