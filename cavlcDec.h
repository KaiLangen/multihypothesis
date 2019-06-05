
#ifndef DECODER_INC_CAVLCDEC_H
#define DECODER_INC_CAVLCDEC_H

#include "defs.h"
#include "cavlc.h"
#include "bitstream.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class CavlcDec : public Cavlc
{
public:
  CavlcDec(Codec* codec, int blockSize);

  int decode(int* iDCT, int ix, int iy, Bitstream* bs);

  void clearNnz(int index) { _mbs[index].nnz[0][0] = 0; };

private:
  int decodeNumCoeffTrailingOnes(int& nc, int& t1s, int vlc, Bitstream* bs);
  int decodeLevel(int iNumCoef,int iTrailingOnes,int *iLevel, Bitstream* bs);
  int decodeTotalZero(int &iTotalZeros,int iNumCoef, Bitstream* bs);
  int decodeRun(int &iRun,int iZerosLeft, Bitstream* bs);
};

#endif // DECODER_INC_CAVLCDEC_H

