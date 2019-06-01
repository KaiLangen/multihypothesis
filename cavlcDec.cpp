
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>

#include "cavlcDec.h"
#include "fileManager.h"
#include "codec.h"
#include "bitstream.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
CavlcDec::CavlcDec(Codec* codec, int blockSize) : Cavlc(codec, blockSize)
{
  for (int i = 0; i < _codec->getBitPlaneLength(); i++)
    _mbs[i].nnz[0][0] = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decode(int* iDCT, int ix, int iy, Bitstream* bs)
{
  int width, numCoeff,t1s,totalZeros,zerosLeft;
  if (bs == _codec->getBitstream())
    width = _codec->getFrameWidth();
  else
    width = _codec->getFrameWidth() >> 1;

  int qp = _codec->getQp();
  int iBitCount=0;
  int nc,vlc;

  int index = (ix/4)+(iy/4)*(width/4);

  numCoeff = t1s = totalZeros = zerosLeft = 0;

  if (ix == 0 && iy == 0)
    nc = 0;
  else if (ix == 0 && iy != 0)
    nc = getNumNonzero(ix, iy-4, width);
  else if (ix != 0 && iy == 0)
    nc = getNumNonzero(ix-4, iy, width);
  else
    nc = (getNumNonzero(ix, iy-4, width) +
          getNumNonzero(ix-4, iy, width) + 1) >> 1;

  if (nc < 2)
    vlc = 0;
  else if (nc < 4)
    vlc = 1;
  else if (nc < 8)
    vlc = 2;
  else
    vlc = 3;

  iBitCount += decodeNumCoeffTrailingOnes(numCoeff, t1s, vlc, bs);

  _mbs[index].nnz[0][0] = numCoeff;

  if (numCoeff == 0)
    return iBitCount;

  iBitCount += decodeLevel(numCoeff, t1s, _mbs[index].level, bs);

  _mbs[index].nnz[0][0] = numCoeff;

  if (numCoeff == 16) {
    totalZeros = 0;
    zerosLeft  = 0;
  }
  else {
    iBitCount += decodeTotalZero(totalZeros, numCoeff, bs);
    zerosLeft  = totalZeros;
  }

  for (int i = 0; i < (numCoeff-1); i++) {
    if (zerosLeft > 0)
      iBitCount += decodeRun(_mbs[index].iRun[i], zerosLeft, bs);
    else
      _mbs[index].iRun[i] = 0;

    zerosLeft -= _mbs[index].iRun[i];
  }

  _mbs[index].iRun[numCoeff-1] = zerosLeft;

  int iSign;
  int iCoeffNum = -1;

#if MODE_DECISION
  iCoeffNum += _codec->getNumChnCodeBands();
#endif

  for (int i = numCoeff-1; i >= 0; i--) {
    iCoeffNum += _mbs[index].iRun[i] + 1;
    assert(iCoeffNum < 16);

    int x = ScanOrder[iCoeffNum][0];
    int y = ScanOrder[iCoeffNum][1];

    iSign = (_mbs[index].level[i] >= 0) ? 0 : 1;
    iDCT[ix+x + (iy+y)*width] = abs(_mbs[index].level[i]);
    if (iSign == 1)
      iDCT[ix+x + (iy+y)*width] |= (0x1 << (_codec->getQuantMatrix(qp, x, y)-1));
  }

  return iBitCount;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeNumCoeffTrailingOnes(int& numCoeff, int& t1s,
                                         int vlc, Bitstream* bs)
{
  int value;

  int length = 0;

  if (vlc < 3) {
    value = 0;

    for (length = 1; length < 17; length++) {
      value <<= 1;
      value |= bs->read(1);

      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 17; i++) {
          if (NumVlcTableC[vlc][j][i] == value &&
              NumVlcTableL[vlc][j][i] == length) {
            numCoeff = i;
            t1s = j;
            goto DecodeNumTrailDone;
          }
        }
    }

    DecodeNumTrailDone:
      ;
  }
  else {
    length = 6;
    value = bs->read(6);

    if (value == 3)
      numCoeff = t1s = 0;
    else {
      numCoeff = (value>>2) + 1;
      t1s = value & 0x3;
    }
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeLevel(int iNumCoef, int iTrailingOnes,
                          int* iLevel, Bitstream* bs)
{
  int iLevelPrefix,iSuffixLength,iLevelSuffixSize;
  unsigned int iLevelSuffix;
  int iLevelCode;
  int iIndex = 0;

  int length = 0;

  //decode trailingOnes
  for (int i = 0; i < iTrailingOnes; i++) {
    int b = bs->read(1);
    length++;

    if (b == 0)
      iLevel[iIndex] = 1;
    else
      iLevel[iIndex] = -1;

    iIndex++;
  }

  //initialize SuffixLength
  if (iNumCoef > 10 && iTrailingOnes < 3)
    iSuffixLength = 1;
  else
    iSuffixLength = 0;

  for(int i = 0; i < iNumCoef-iTrailingOnes; i++) {
    iLevel[iIndex] = 0;

    //get level prefix
    iLevelPrefix = -1;

    for (int b = 0; !b; iLevelPrefix++) {
      b = bs->read(1);
      length++;
    }

    if (iLevelPrefix == 14 && iSuffixLength == 0)
      iLevelSuffixSize = 4;
    else if (iLevelPrefix >= 15)
      iLevelSuffixSize = iLevelPrefix - 3;
    else
      iLevelSuffixSize = iSuffixLength;

    if (iLevelSuffixSize > 0) {
      iLevelSuffix = (unsigned int)bs->read(iLevelSuffixSize);
      length += iLevelSuffixSize;
    }
    else
      iLevelSuffix = 0;

    iLevelCode = (Min(15,iLevelPrefix)<<iSuffixLength) + iLevelSuffix;

    if (iLevelPrefix >= 15 && iSuffixLength == 0)
      iLevelCode += 15;
    if (iLevelPrefix >= 16)
      iLevelCode += ((1<<( iLevelPrefix-3 ))-4096);
    if (iIndex == iTrailingOnes && iTrailingOnes < 3)
      iLevelCode += 2;

    if (iLevelCode % 2 == 0) //even number
      iLevel[iIndex] = (iLevelCode + 2 ) >> 1;
    else
      iLevel[iIndex] = (-iLevelCode-1) >> 1;

    if (iSuffixLength == 0)
      iSuffixLength = 1;
    if((abs(iLevel[iIndex]) > (3<<(iSuffixLength-1))) && (iSuffixLength < 6))
      iSuffixLength++;

    iIndex++;
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeTotalZero(int& iTotalZeros, int iNumCoef, Bitstream* bs)
{
  int iValue = 0;

  int length = 0;

  for (length = 1; length < 10; length++) {
    int success = 0;
    iValue <<= 1;
    iValue |= bs->read(1);

    for (int i = 0; i < 16; i++) {
      if (TotalZerosTableC[iNumCoef-1][i] == iValue &&
          TotalZerosTableL[iNumCoef-1][i] == length) {
        success = 1;
        iTotalZeros = i;
        break;
      }
    }

    if (success == 1)
      break;
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeRun(int& iRun, int iZerosLeft, Bitstream* bs)
{
  int iValue = 0;

  int length = 0;

  for (length = 1; length < 10; length++) {
    int success = 0;
    iValue <<= 1;
    iValue |= bs->read(1);

    for (int i = 0; i < 15; i++) {
      if (RunTableC[iZerosLeft-1][i] == iValue &&
          RunTableL[iZerosLeft-1][i] == length) {
        success = 1;
        iRun = i;
        break;
      }
    }

    if (success == 1)
      break;
  }

  return length;
}

