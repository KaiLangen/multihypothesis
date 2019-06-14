
#ifndef DECODER_INC_SIDEINFORMATION_H
#define DECODER_INC_SIDEINFORMATION_H

#include <map>

#include "defs.h"
#include "frameBuffer.h"

class CorrModel;
class Codec;

class SideInformation
{
public:
  SideInformation(Codec* codec, CorrModel* model,
                  std::map<std::string, std::string> configMap);

  void sideInfoMCI(imgpel* prevFrame, imgpel* nextFrame, imgpel* currFrame);

  void chroma_MEMC(RefBuffer* refFrames, imgpel* sideInfo);

  void getResidualFrame(imgpel* bRefFrame, imgpel* fRef,
                        imgpel* curr, int* residue, int* rcList);

  void getRecFrame(imgpel *imgBRef, imgpel *imgFRef,
                   int *iResidue, imgpel *imgRec, int *rcList);

  void getRefinedSideInfo(imgpel *imgPrevKey, imgpel *imgNextKey,
                          imgpel *imgCurrFrame, imgpel* imgTmpRec,
                          imgpel *imgRefined, int iMode);

private:
  void lowpassFilter(imgpel* src, imgpel* dst, const int boxSize);

  // MCI functions
  void createSideInfoProcess(imgpel* prevFrame, imgpel* nextFrame,
                             imgpel* MCfwd, imgpel* MCback, int iMode);

  void forwardME(imgpel* prev, imgpel* curr,
                 mvinfo* candidate, const int iRange);

  void getRefinedSideInfoProcess(imgpel* imgPrevBuffer, imgpel* imgTmpRec,
                                 imgpel* imgSI, imgpel* imgRefined,
                                 mvinfo* varList, int iMode);

  void bidirectME(imgpel* imgPrev, imgpel* imgNext, mvinfo* varCandidate,
                  const int iPadSize, const int iRange);

  void spatialSmooth(imgpel* imgPrev, imgpel* imgNext, mvinfo* varCandidate,
                     const int iBlockSize, const int iPadSize);

  void MC(imgpel* imgPrev, imgpel* imgNext, imgpel* imgDst,
          imgpel* imgMCf, imgpel* imgMCb,
          mvinfo* varCandidate, mvinfo* varCandidate2,
          const int iPadSize, const int iRange, const int iMode);


  // Chroma-ME functions
  void TSS(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
           mvinfo& mv, int step, int center);

  void ES(imgpel* trgU, imgpel* trgV, imgpel* refU, imgpel* refV,
           mvinfo& mv, int p, int center, int padSize);

  void ME(imgpel* refFrameU, imgpel* currFrameU,
          imgpel* refFrameV, imgpel* currFrameV,
          mvinfo* candidate);

  void MC(imgpel* sideInfo, vector<mvinfo*> mvs,
          vector<imgpel*> refs, int padSize);

  void MC(imgpel* sideInfo, mvinfo* candidate,
          imgpel* ref, int padSize);

  void MC(imgpel* refs[], mvinfo* mvs[],
          imgpel* imgDst, int padSize, int nRefs);

  void spatialSmooth(imgpel* rU, imgpel* rV,
                     imgpel* cU, imgpel* cV,
                     mvinfo* varCandidate,
                     const int iBlockSize,
                     const int iPadSize);

  // general
  void getSkippedRecFrame(imgpel* prevKey, imgpel* imgWZFrame, int* skipMask);

  Codec* _codec;
  CorrModel* _model;

  int _width;
  int _height;
  int _frameSize;
  int _blockSize;
  int _ss;
  int _p;
  int _nmv;
  int* _skipMask;

  // SI refinement
  int* _refinedMask;
  mvinfo*     _varList0;
  mvinfo*     _varList1;

# if OBMC
  const static int _H[3][8][8];
# endif
};

float getDCValue(imgpel* img,int iStep,int iStripe,int iBlock);


#endif // DECODER_INC_SIDEINFORMATION_H

