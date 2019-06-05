
#ifndef DECODER_INC_DECODER_H
#define DECODER_INC_DECODER_H

#include <iostream>

#include "defs.h"
#include "codec.h"

using namespace std;

class FileManager;
class Transform;
class CorrModel;
class SideInformation;
class CavlcDec;
class FrameBuffer;
class LdpcaDec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Decoder : public Codec
{
public:
  Decoder(map<string, string> configMap);
  ~Decoder() { /* TODO Remember to free memory space */ };

  void decodeWzFrame();

  int* getSpiralSearchX()     { return _spiralSearchX; };
  int* getSpiralSearchY()     { return _spiralSearchY; };
  int* getSpiralHpelSearchX() { return _spiralHpelSearchX; };
  int* getSpiralHpelSearchY() { return _spiralHpelSearchY; };
  int* getSkipMask()          { return _skipMask; };

  int* _skipMask;

private:
  void initialize();

  void decodeWzHeader();

  void parseKeyStat(const char* filename, double* rate, double* psnr);

  int getSyndromeData();
  int decodeSkipMask();

  void getSkippedRecFrame(imgpel* imgPrevFrame, imgpel* imgWZFrame, int* skipMask);

  void getSourceBit(int* dct_q, double* source, int q_i, int q_j, int curr_pos);
  double decodeLDPC(int* iQuantDCT, int* iDCT, int* iDecoded, int x, int y, int iOffset);

  void motionSearchInit(int maxsearch_range);

private:
  FileManager*      _files;

  FrameBuffer*      _fb;

  Transform*        _trans;

  CorrModel*        _model;
  SideInformation*  _si;

  CavlcDec*         _cavlc;
  CavlcDec*         _cavlcU;
  CavlcDec*         _cavlcV;
  LdpcaDec*         _ldpca;

  int               _maxValue[4][4];

  int               _rcBitPlaneNum;
  int               _rcQuantMatrix[4][4];

  int*              _spiralSearchX;
  int*              _spiralSearchY;
  int*              _spiralHpelSearchX;
  int*              _spiralHpelSearchY;
};

void decodeBits(double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors,unsigned char crccode,int numcode);
int  beliefPropagation(int *ir, int *jc, int m, int n, int nzmax,
                       double *LLR_intrinsic, double *syndrome,
                       double *decoded);

bool checkCRC(double * source,const int length,unsigned char crc);

int getSymbol(int len,int &curr_pos,char *buffer);

#endif // DECODER_INC_DECODER_H

