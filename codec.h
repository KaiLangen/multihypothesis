
#ifndef COMMON_INC_CODEC_H
#define COMMON_INC_CODEC_H

#include <map>

#include "defs.h"

class Bitstream;

class Codec
{
public:
  Codec() {};

  int     getFrameWidth()      { return _frameWidth; };
  int     getFrameHeight()     { return _frameHeight; };
  int     getBitPlaneLength()  { return _bitPlaneLength; };
  int     getQp()              { return _qp; };
  int     getChrQp()              { return _chrQp; };
  int     getKeyQp()           { return _keyQp; };
  int     getNumChnCodeBands() { return _numChnCodeBands; };

  double* getAverage()         { return _average; };
  double* getAlpha()           { return _alpha; };
  double* getSigma()           { return _sigma; };

  int     getQuantMatrix(int qp, int x, int y){return QuantMatrix[qp][y][x];};
  int     getQuantStep(int x, int y, bool isChr)
  { if (isChr) return _qStepChr[y][x];
          else return _quantStep[y][x];};

  void getRecFrame(imgpel* recon, imgpel* bRef, imgpel* fRef,
                   int* curr,  int* rcList, bool isChr);
  Bitstream* getBitstream() { return _bs; };

protected:
  const static int  ResidualBlockSize;
  const static int  SkipBlockSize;

  const static int  QuantMatrix[8][4][4];
  const static int  BitPlaneNum[8];
  const static int  MaxBitPlane[4][4];
  const static int  MinQStepSize[8][4][4];

  const static int  ScanOrder[16][2];
  const static int  HuffmanCodeValue[4][3][16];
  const static int  HuffmanCodeLength[4][3][16];

  int               _quantStep[4][4];
  int               _qStepChr[4][4];

  int               _frameWidth;
  int               _frameHeight;
  int               _frameSize;
  int               _numFrames;
  int               _bitPlaneLength;
  int               _qp;
  int               _chrQp;
  int               _keyQp;
  int               _gop;
  int               _numChnCodeBands;

  bool*             _parity;
  double*           _dParity; // TODO temporary for decoder
  unsigned char*    _crc;
  double*           _average;
  double*           _alpha;
  double*           _sigma;

  Bitstream*        _bs;
  Bitstream*        _bsU;
  Bitstream*        _bsV;
};

std::map<std::string,std::string>& readConfig(std::string filename,bool isEnc);

#endif // COMMON_INC_CODEC_H

