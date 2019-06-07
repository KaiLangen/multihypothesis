
#ifndef DECODER_INC_CORRMODEL_H
#define DECODER_INC_CORRMODEL_H

#include "defs.h"

class Transform;
class Codec;

class CorrModel
{
public:
  CorrModel(Codec* codec, Transform* trans) { _codec = codec; _trans = trans; };

  void updateCNM(imgpel* imgForward,imgpel* imgBackward,int *g_iRefinedMask);

  void correlationNoiseModeling(imgpel* imgMCForward,imgpel *imgMCBackward);

#if SKIP_MODE
double getSoftInput(int* si, int* skipMask, int iCurrPos, int *iDecoded,
                    double* dLLR, int x, int y, int iMode);
#else
double getSoftInput(int* si, int iCurrPos, int *iDecoded,
                    double* dLLR, int x, int y, int iMode);
#endif

private:
  Codec*      _codec;

  Transform*  _trans;
};

double getLaplacian(int iLowerBound,int iUpperBound,double dAlpha,int iSideInfo,int iQuantSize,int iMode);
double getCondProb(int iBit,int iBand,int iBitLength,int iCurrPos,int iDecodedBit,double dAlpha,int iQuantSize,int iMode);

#endif // DECODER_INC_CORRMODEL_H

