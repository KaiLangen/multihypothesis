#include "frameBuffer.h"
#include "sideInformation.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* RefBuffer
*/
RefBuffer::RefBuffer(int width, int height, int windowSize)
{
  _windowSize = windowSize;
  _nextRec = 0;
  _isFull = false;
  _padSize = 40;
  _width = width;
  _height = height;
  _frameSize = width * height;
  _yuvFrameSize = 3*_frameSize>>1;
  _paddedFrameSize = (width+2*_padSize) * (height+2*_padSize);
  _refFrames.resize(windowSize, vector<imgpel*>(4));
  for (size_t i = 0; i < _refFrames.size(); i++) {
    _refFrames[i][0] = new imgpel[_yuvFrameSize];
    _refFrames[i][1] = new imgpel[_paddedFrameSize];
    _refFrames[i][2] = new imgpel[_paddedFrameSize];
    _refFrames[i][3] = new imgpel[_paddedFrameSize];
  }
  _currFrame.push_back(new imgpel[_yuvFrameSize]);
  _currFrame.push_back(new imgpel[_paddedFrameSize]);
  _currFrame.push_back(new imgpel[_paddedFrameSize]);
  _currFrame.push_back(new imgpel[_paddedFrameSize]);
}

RefBuffer::~RefBuffer()
{
  for (size_t v = 0; v < _refFrames.size(); v++) {
    for (int i = 0; i < 4; i++)
      delete [] _refFrames[v][i];
  }

  // delete curr
  for (int i = 0; i < 4; i++) {
    delete [] _currFrame[i];
  }
}

void RefBuffer::updateRecWindow(FILE* fKeyReadPtr, int prevKeyNo)
{
  // loop to start of circular buffer
  if (_nextRec >= _windowSize) {
    _nextRec = 0;
    _isFull = true;
  }
  // read the newest frame from file
  imgpel* newKey = _refFrames[_nextRec][0];
  fseek(fKeyReadPtr, (prevKeyNo+_windowSize)*_yuvFrameSize, SEEK_SET);
  fread(newKey, _yuvFrameSize, 1, fKeyReadPtr);
  int ww = _width>>1;
  int hh = _height>>1;
  imgpel* UChroma = new imgpel[_frameSize];
  imgpel* VChroma = new imgpel[_frameSize];

  // upsample Chroma planes and pad
  bilinear(newKey+_frameSize, UChroma,
           ww, hh, ww, hh, 0, 0);
  bilinear(newKey+5*(_frameSize>>2), VChroma,
           ww, hh, ww, hh, 0, 0);
  
  pad(newKey, _refFrames[_nextRec][3], _width, _height, 40);
  pad(UChroma, _refFrames[_nextRec][1], _width, _height, 40);
  pad(VChroma, _refFrames[_nextRec][2], _width, _height, 40);

  _nextRec++;
  delete [] UChroma;
  delete [] VChroma;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void RefBuffer::init(FILE* fKeyReadPtr) {
  int ww = _width>>1;
  int hh = _height>>1;
  int chsize = _frameSize>>2;

  imgpel* UChroma = new imgpel[_frameSize];
  imgpel* VChroma = new imgpel[_frameSize];

  _refFrames.resize(_windowSize, vector<imgpel*>(4));
  // allocate space for reference frames
  for (size_t i = 0; i < _refFrames.size(); i++) {
    _refFrames[i][0] = new imgpel[_yuvFrameSize];
    _refFrames[i][1] = new imgpel[_paddedFrameSize];
    _refFrames[i][2] = new imgpel[_paddedFrameSize];
    _refFrames[i][3] = new imgpel[_paddedFrameSize];
  }

  // read data from file and pad into extra buffers
  for (size_t i = 0; i < _refFrames.size(); i++) {
    fseek(fKeyReadPtr, _yuvFrameSize*i, SEEK_SET);
    fread(_refFrames[i][0], _yuvFrameSize, 1, fKeyReadPtr);
    imgpel* ptrU = _refFrames[i][0] + _frameSize;
    imgpel* ptrV = ptrU + chsize;

    bilinear(ptrU, UChroma, ww, hh, ww, hh, 0, 0);
    bilinear(ptrV, VChroma, ww, hh, ww, hh, 0, 0);
    pad(UChroma, _refFrames[i][1], _width, _height, 40);
    pad(VChroma, _refFrames[i][2], _width, _height, 40);
    pad(_refFrames[i][0], _refFrames[i][3], _width, _height, 40);
    _nextRec++;
  }

  delete [] UChroma;
  delete [] VChroma;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* FrameBuffer
*/
FrameBuffer::FrameBuffer(int width, int height, int windowSize)
{
  _frameSize        = width * height;
  _origFrame        = new imgpel[3*(_frameSize>>1)];
  _sideInfoFrame    = new imgpel[3*(_frameSize>>1)];
  _rBuff            = new RefBuffer(width, height, windowSize);

}

FrameBuffer::~FrameBuffer()
{
  delete [] _origFrame;
  delete [] _sideInfoFrame;
  delete _rBuff;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* General frame modification functions
*/
void pad(imgpel* src, imgpel* dst, int width, int height, int iPadSize)
{
  int padWidth  = width  + 2*iPadSize;
  int padHeight = height + 2*iPadSize;

  // Loops start from iPadSize; subtract it from src index

  // Upper left
  for (int y = 0; y < iPadSize; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padWidth] = src[0];

  // Upper
  for (int x = iPadSize; x < iPadSize+width; x++)
    for (int y = 0; y < iPadSize; y++)
      dst[x+y*padWidth] = src[x-iPadSize];

  // Upper right
  for (int y = 0; y < iPadSize; y++)
    for (int x = iPadSize+width; x < padWidth; x++)
      dst[x+y*padWidth] = src[width-1];

  // Left
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padWidth] = src[(y-iPadSize)*width];

  // Middle
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = iPadSize; x < iPadSize+width; x++)
      dst[x+y*padWidth] = src[x-iPadSize+(y-iPadSize)*width];

  // Right
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = iPadSize+width; x < padWidth; x++)
      dst[x+y*padWidth] = src[(y-iPadSize+1)*width-1];

  // Bottom left
  for (int y = iPadSize+height; y < padHeight; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*padWidth] = src[(height-1)*width];

  // Bottom
  for (int y = iPadSize+height; y < padHeight; y++)
    for (int x = iPadSize; x < iPadSize+width; x++)
      dst[x+y*padWidth] = src[(x-iPadSize)+(height-1)*width];

  // Bottom right
  for (int y = iPadSize+height; y < padHeight; y++)
    for (int x = iPadSize+width; x < padWidth; x++)
      dst[x+y*padWidth] = src[height*width-1];
}


void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight, int px, int py){
  for(int j=0;j<buffer_h;j++)
    for(int i=0;i<buffer_w;i++)
    {
      int buffer_r=2*buffer_w;
      int a,b,c,d;

      int x=px+i;
      int y=py+j;
      if(x>picwidth-1)x=picwidth-1;
      if(x<0)x=0;
      if(y>picheight-1)y=picheight-1;
      if(y<0)y=0;
      a=source[(x)+(y)*picwidth];

      if((x+1)<picwidth) b=source[(x+1)+(y)*picwidth];
      else b=a;

      if((y+1)<picheight) c=source[(x)+(y+1)*picwidth];
      else c=a;

      if((x+1)<picwidth && (y+1)<picheight) d=source[(x+1)+(y+1)*picwidth];
      else d=a;

      buffer[2*i+(2*j)*buffer_r]=a;
      buffer[(2*i+1)+(2*j)*buffer_r]=(a+b)/2;
      buffer[2*i+(2*j+1)*buffer_r]=(a+c)/2;
      buffer[(2*i+1)+(2*j+1)*buffer_r]=(a+b+c+d)/4;
    }
}

void decimate2(imgpel* source, imgpel* buffer, int buffer_w, int buffer_h,
               int picwidth, int px, int py){
  for(int j=0;j<buffer_h;j++)
    for(int i=0;i<buffer_w;i++)
      buffer[i + j*buffer_w] = source[(px+2*i) + (py+2*j)*picwidth];
}
