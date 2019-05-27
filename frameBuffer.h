
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include "defs.h"

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int gop = 0)
  {
    (void)gop;
    int frameSize = width * height;
    _prevFrame        = new imgpel[frameSize];
    _currFrame        = new imgpel[frameSize];
    _nextFrame        = new imgpel[frameSize];
    _prevChroma       = new imgpel[frameSize>>1];
    _currChroma       = new imgpel[frameSize>>1];
    _sideInfoFrame    = new imgpel[frameSize];
    _origFrame        = new imgpel[frameSize];

//    if (gop != 0) {
//      _refFrames      = new imgpel*[gop-1];
//
//      for (int i = 0; i < gop-1; i++)
//        _refFrames[i] = new imgpel[frameSize];
//    }

    _dctFrame         = new int[frameSize];
    _quantDctFrame    = new int[frameSize];
    _decFrame         = new int[frameSize];
    _invQuantDecFrame = new int[frameSize];
  };

  imgpel*  getPrevFrame()        { return _prevFrame; };
  imgpel*  getPrevChroma()       { return _prevChroma; };
  imgpel*  getCurrFrame()        { return _currFrame; };
  imgpel*  getCurrChroma()       { return _currChroma; };
  imgpel*  getNextFrame()        { return _nextFrame; };
  imgpel*  getorigFrame()        { return _origFrame; };
  imgpel*  getSideInfoFrame()    { return _sideInfoFrame; };
  int*     getDctFrame()         { return _dctFrame; };
  int*     getQuantDctFrame()    { return _quantDctFrame; };
  int*     getDecFrame()         { return _decFrame; };
  int*     getInvQuantDecFrame() { return _invQuantDecFrame; };

private:
  imgpel*  _prevFrame;
  imgpel*  _prevChroma;
  imgpel*  _currFrame;
  imgpel*  _currChroma;
  imgpel*  _nextFrame;
  imgpel*  _origFrame;
  imgpel*  _sideInfoFrame;
  int*     _dctFrame;
  int*     _quantDctFrame;
  int*     _decFrame;
  int*     _invQuantDecFrame;
};

#endif // ENCODER_INC_FRAMEBUFFER_H

