
#ifndef ENCODER_INC_ENCODER_H
#define ENCODER_INC_ENCODER_H

#include "defs.h"
#include "codec.h"

class FileManager;
class FrameBuffer;
class Transform;
class CavlcEnc;
class LdpcaEnc;
class Codec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Encoder : public Codec
{
public:
  Encoder(std::map<std::string, std::string> configMap);
  ~Encoder() { /* TODO Remember to free memory space */ };

  void encodeKeyFrame();
  void encodeWzFrame();
  void encodeWzHeader();

protected:
  void initialize();

  void updateMaxValue(int* block);

  void selectCodingMode(int* frame);

  void generateSkipMask();

  int encodeSkipMask();
  int getHuffmanCode(int qp, int type, int symbol, int& code, int& length);

  void encodeFrameLdpca(int* frame);
  void setupLdpcaSource(int* frame, int* source, int offsetX, int offsetY, int bitPosition);
  void computeCRC(int* data, const int length, unsigned char* crc);

  void report();

protected:
  const static int  Scale[3][8];

  FileManager*      _files;

  FrameBuffer*      _fb;

  Transform*        _trans;

  CavlcEnc*         _cavlc;
  CavlcEnc*         _cavlcU;
  CavlcEnc*         _cavlcV;
  LdpcaEnc*         _ldpca;

  int               _rcBitPlaneNum;
  int               _rcQuantMatrix[4][4];

  int               _maxValue[4][4];
  int*              _skipMask;
  int               _prevMode;
  int               _prevType;

  int               _modeCounter[4];
  unsigned char*    _crcPtr;
};

#endif // ENCODER_INC_ENCODER_H

