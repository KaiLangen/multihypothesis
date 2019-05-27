
#ifndef ENCODER_INC_CAVLCENC_H
#define ENCODER_INC_CAVLCENC_H

#include <vector>

#include <cstdio>

#include "defs.h"
#include "cavlc.h"
#include "bitstream.h"

using std::vector;

class File;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class SyntaxElement
{
public:
  SyntaxElement(int value1 = 0, int value2 = 0, int length = 0, int info = 0) {
    _value1 = value1;
    _value2 = value2;
    _length = length;
    _info   = info;
  };

  void set(int value1, int value2, int length, int info) {
    _value1 = value1;
    _value2 = value2;
    _length = length;
    _info   = info;
  };

  int getValue1() { return _value1; };
  int getValue2() { return _value2; };
  int getLength() { return _length; };
  int getInfo()   { return _info; };

  void setValue1(int value1) { _value1 = value1; };
  void setValue2(int value2) { _value2 = value2; };
  void setLength(int length) { _length = length; };
  void setInfo  (int info)   { _info   = info; };

private:
  int _value1;  //!< numerical value of syntax element
  int _value2;  //!< for blocked symbols, e.g. run/level
  int _length;  //!< length of code
  int _info;    //!< info part of UVLC code
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class CavlcEnc : public Cavlc
{
public:
  CavlcEnc(Codec* codec, int blockSize);

  int encode(int* frame, int* skipMask);
  int encode(int* frame, Bitstream* bs);

private:
  void setupMacroBlock(int* frame, int mbX, int mbY, int width);
  int encodeMacroBlock(int mbX, int mbY, int width, Bitstream* bs);

  int symbol2vlc(SyntaxElement* sym, Bitstream* bs);

  int encodeNumTrail(SyntaxElement* se, Bitstream* bs);
  int encodeSignTrail(vector<int>& sign, Bitstream* bs);
  int encodeLevelsVlc0(SyntaxElement* se, Bitstream* bs);
  int encodeLevelsVlcN(SyntaxElement* se, int vlc, Bitstream* bs);
  int encodeTotalZeros(SyntaxElement* se, Bitstream* bs);
  int encodeRuns(SyntaxElement* se, Bitstream* bs);

};

#endif // ENCODER_INC_CAVLCENC_H

