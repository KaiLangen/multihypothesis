
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
  void getResidualFrame(imgpel* bRefFrame, imgpel* currFrame, int* residue);
  void getRecFrame(imgpel *imgBReference, int *iResidue, imgpel *imgRec);

private:
  void lowpassFilter(imgpel* src, imgpel* dst, const int boxSize);

  int TSS(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
          mvinfo& mv, int step, int center);

  void ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
          mvinfo& mv, int p, int center);

  void ME(imgpel* refFrameU, imgpel* currFrameU,
          imgpel* refFrameV, imgpel* currFrameV);

  void MC(imgpel* imgPrev, imgpel* imgDst);

  void spatialSmooth(imgpel* rU, imgpel* rV, imgpel* cU, imgpel* cV, mvinfo* varCandidate,
                     const int iBlockSize, const int iPadSize);

  void pad(imgpel* src, imgpel* dst, const int iPadSize);

  bool isSkip(const int* skipMask, int start, int thresh);

  Codec*      _codec;
  CorrModel*  _model;

  int _width;
  int _height;
  int _frameSize;
  int _blockSize;
  int _p;
  int _nmv;

# if OBMC
  const static int _H[3][8][8];
# endif

  mvinfo* _mvs;
  int* _skipMask;
};


void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight);

int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2,
            int s1,int s2,int blocksize);
#endif // DECODER_INC_SIDEINFORMATION_H

