
#ifndef DECODER_INC_SIDEINFORMATION_H
#define DECODER_INC_SIDEINFORMATION_H

#include "defs.h"

class CorrModel;
class Codec;

class SideInformation
{
public:
  SideInformation(Codec* codec, CorrModel* model);

  void createSideInfo(imgpel* prevChroma, imgpel* currChroma,
                      imgpel* imgPrevKey, imgpel* imgCurrFrame);
# if RESIDUAL_CODING
  void getResidualFrame(imgpel* bRefFrame, imgpel* currFrame, int* residue);
  void getRecFrame(imgpel *imgBReference, int *iResidue, imgpel *imgRec);
# endif

private:
  void lowpassFilter(imgpel* src, imgpel* dst, const int boxSize);

  void ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
          mvinfo& mvs, int p, int center);

  void ME(imgpel* refFrameU, imgpel* currFrameU,
          imgpel* refFrameV, imgpel* currFrameV,
          mvinfo* mvs);

  void MC(imgpel* imgPrev, imgpel* imgDst);

  Codec*      _codec;
  CorrModel*  _model;

  int _width;
  int _height;
  int _frameSize;
  int _blockSize;
  int _nmv;

# if OBMC
  const static int _H[3][8][8];
  void OBMC(imgpel* imgPrev, imgpel* imgDst);
# endif

  mvinfo* _mvs;
};

#endif // DECODER_INC_SIDEINFORMATION_H

