
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include <vector>
#include <cstring>

#include "defs.h"

using namespace std;
class RefBuffer
{
public:
  vector<vector<imgpel*>> _refFrames;
  vector<imgpel*> _currFrame;
  vector<imgpel*> _prevKeyFrame;
  vector<imgpel*> _nextKeyFrame;
  int _windowSize;
  int _frameSize;
  int _yuvFrameSize;
  int _paddedSize;
  int _nextRec;
  bool _isFull;

  RefBuffer(int width, int height, int windowSize)
  {
    _windowSize = windowSize;
    _nextRec = 0;
    _isFull = false;
    _frameSize = width * height;
    _yuvFrameSize = 3*_frameSize>>1;
    _paddedSize = (width+80) * (height+80);
    _refFrames.resize(windowSize, vector<imgpel*>(4));
    for (auto v : _refFrames) {
      v[0] = new imgpel[_yuvFrameSize];
      v[1] = new imgpel[_paddedSize];
      v[2] = new imgpel[_paddedSize];
      v[3] = new imgpel[_paddedSize];
    }
    _currFrame.push_back(new imgpel[_yuvFrameSize]);
    _currFrame.push_back(new imgpel[_paddedSize]);
    _currFrame.push_back(new imgpel[_paddedSize]);
    _currFrame.push_back(new imgpel[_paddedSize]);
    _prevKeyFrame.push_back(new imgpel[_yuvFrameSize]);
    _prevKeyFrame.push_back(new imgpel[_paddedSize]);
    _prevKeyFrame.push_back(new imgpel[_paddedSize]);
    _prevKeyFrame.push_back(new imgpel[_paddedSize]);
    _nextKeyFrame.push_back(new imgpel[_yuvFrameSize]);
    _nextKeyFrame.push_back(new imgpel[_paddedSize]);
    _nextKeyFrame.push_back(new imgpel[_paddedSize]);
    _nextKeyFrame.push_back(new imgpel[_paddedSize]);
  }

  ~RefBuffer()
  {
    for (auto v : _refFrames)
      for (int i = 0; i < 4; i++)
        delete [] v[i];

    // delete curr, prev, next
    for (int i = 0; i < 4; i++) {
      delete [] _currFrame[i];
      delete [] _prevKeyFrame[i];
      delete [] _nextKeyFrame[i];
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  vector<vector<imgpel*>>::iterator begin() { return _refFrames.begin(); }

  vector<vector<imgpel*>>::iterator end()
  {
    if (_isFull) return _refFrames.end();
    else return _refFrames.begin() + _nextRec;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  vector<imgpel*> getCurrFrame() { return _currFrame; }


  // ---------------------------------------------------------------------------
  // increment nextRec, and wrap if equal to buffSize
  // ---------------------------------------------------------------------------
  void updateRecWindow()
  {
    if (_windowSize == 0) return;
    if (_nextRec >= _windowSize) {
      _nextRec = 0;
      _isFull = true;
    }
    memcpy(_refFrames[_nextRec][0], _currFrame[0], _yuvFrameSize);
    memcpy(_refFrames[_nextRec][1], _currFrame[1], _paddedSize);
    memcpy(_refFrames[_nextRec][2], _currFrame[2], _paddedSize);
    memcpy(_refFrames[_nextRec][3], _currFrame[3], _paddedSize);
    _nextRec++;
  }

};

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int windowSize = 0)
  {
    _frameSize        = width * height;
    _origFrame        = new imgpel[3*(_frameSize>>1)];
    _sideInfoFrame    = new imgpel[3*(_frameSize>>1)];
    _dctFrame         = new int[3*(_frameSize>>1)];
    _quantDctFrame    = new int[3*(_frameSize>>1)];
    _decFrame         = new int[3*(_frameSize>>1)];
    _invQuantDecFrame = new int[3*(_frameSize>>1)];
    _rBuff            = new RefBuffer(width, height, windowSize);

  };

  ~FrameBuffer()
  {
    delete [] _origFrame;
    delete [] _sideInfoFrame;
    delete [] _dctFrame;
    delete [] _quantDctFrame;
    delete [] _decFrame;
    delete [] _invQuantDecFrame;
    delete _rBuff;
  }

  //  Luma
  imgpel*  getPrevFrame()        { return _rBuff->_prevKeyFrame[0]; };
  imgpel*  getCurrFrame()        { return _rBuff->_currFrame[0]; };
  imgpel*  getNextFrame()        { return _rBuff->_nextKeyFrame[0]; };
  imgpel*  getorigFrame()        { return _origFrame; };
  imgpel*  getSideInfoFrame()    { return _sideInfoFrame; };
  int*     getDctFrame()         { return _dctFrame; };
  int*     getQuantDctFrame()    { return _quantDctFrame; };
  int*     getDecFrame()         { return _decFrame; };
  int*     getInvQuantDecFrame() { return _invQuantDecFrame; };
  RefBuffer* getRefBuffer()      { return _rBuff; };

private:
  int      _frameSize;
  imgpel*  _origFrame;
  imgpel*  _sideInfoFrame;

  int*     _dctFrame;
  int*     _quantDctFrame;
  int*     _decFrame;
  int*     _invQuantDecFrame;
  RefBuffer* _rBuff;
};

#endif // ENCODER_INC_FRAMEBUFFER_H

