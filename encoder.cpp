
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "config.h"
#include "encoder.h"
#include "transform.h"
#include "fileManager.h"
#include "frameBuffer.h"
#include "cavlcEnc.h"
#include "ldpcaEnc.h"
#include "bitstream.h"

using namespace std;

const int Encoder::Scale[3][8] = {
  {8, 6, 6, 4, 4, 3, 2, 1},
  {8, 8, 8, 4, 4, 4, 2, 1},
  {4, 4, 4, 4, 3, 2, 2, 1}
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Encoder::Encoder(map<string, string> configMap)
{
  _files = FileManager::getManager();

  _qp          = atoi(configMap["WzQP"].c_str()); 
  _keyQp       = atoi(configMap["KeyQP"].c_str());
  if (!strcmp(configMap["SequenceType"].c_str(), "CIF")) {
    _frameWidth  = 352;
    _frameHeight = 288;
  } else {
    _frameWidth  = 176;
    _frameHeight = 144;
  }
  _numFrames     = atoi(configMap["NumFrames"].c_str());
  _gop           = atoi(configMap["Gop"].c_str());

  string wzFileName = configMap["WZFile"];
  _files->addFile("src", configMap["SrcFile"])->openFile("rb");
  _files->addFile("wz",  wzFileName)->openFile("wb");
  _files->addFile("key", configMap["KeyFile"]);
  _files->addFile("chroma", configMap["ChromaFile"]);

  string ubs = wzFileName.substr(0, wzFileName.find(".bin")) + ".u.bin";
  string vbs = wzFileName.substr(0, wzFileName.find(".bin")) + ".v.bin";
  _files->addFile("wzU",  ubs.c_str())->openFile("wb");
  _files->addFile("wzV",  vbs.c_str())->openFile("wb");

  _bs = new Bitstream(1024, _files->getFile("wz")->getFileHandle());
  _bsU = new Bitstream(1024, _files->getFile("wzU")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV")->getFileHandle());
  initialize();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::initialize()
{
  _frameSize        = _frameWidth * _frameHeight;
  _bitPlaneLength   = _frameSize / 16;

  _numChnCodeBands  = 16;

  _parity           = new bool[_bitPlaneLength * BitPlaneNum[_qp]];

# if HARDWARE_LDPC
  if (_bitPlaneLength == 6336)
    _crc            = new unsigned char[BitPlaneNum[_qp] * 4];
  else
    _crc            = new unsigned char[BitPlaneNum[_qp]];
# else
  _crc              = new unsigned char[BitPlaneNum[_qp]];
# endif
  _crcPtr           = _crc;

  _average          = new double[16];
  _alpha            = new double[_frameSize];
  _sigma            = new double[16];

  _skipMask         = new int[_bitPlaneLength];

  _prevMode         = 0;
  _prevType         = 0;

  for (int i = 0; i < 4; i++)
    _modeCounter[i] = 0;

  _fb = new FrameBuffer(_frameWidth, _frameHeight);

  _trans = new Transform(this);

  _cavlc = new CavlcEnc(this, 4);
  _cavlcU = new CavlcEnc(this, 4);
  _cavlcV = new CavlcEnc(this, 4);

  // Initialize LDPC
  string ladderFile;

# if HARDWARE_LDPC
  ladderFile = "ldpca/1584_regDeg3.lad";
# else
  if (_frameWidth == 352 && _frameHeight == 288)
    ladderFile = "ldpca/6336_regDeg3.lad";
  else
    ladderFile = "ldpca/1584_regDeg3.lad";
# endif

  _ldpca = new LdpcaEnc(ladderFile, this);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeKeyFrame()
{
  string srcFileName = _files->getFile("src")->getFileName();
  string keyFileName = _files->getFile("key")->getFileName();

  cout << "Running JM to encode key frames" << endl;

  stringstream cmd(stringstream::in | stringstream::out);
  cmd << "cd jm; ";
  cmd << "./lencod.exe -d encoder_intra_main.cfg ";
  cmd << "-p InputFile=\"" << BIN_DIR << "/" << srcFileName << "\" ";
  cmd << "-p ReconFile=\"" << BIN_DIR << "/" << keyFileName << "\" ";
  cmd << "-p FramesToBeEncoded=" << ((_numFrames + _gop/2)/_gop) << " ";
  cmd << "-p QPISlice=" << _keyQp << " ";
  cmd << "-p FrameSkip=" << _gop-1 << " ";
  cmd << "-p SourceWidth=" << _frameWidth << " ";
  cmd << "-p SourceHeight=" << _frameHeight << " ";
  cmd << " > jm.log;";
  cmd << "cp stats.dat " << BIN_DIR;
  system(cmd.str().c_str());
  _files->getFile("key")->openFile("rb");
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeWzHeader()
{
  // Encode width and height in terms of macroblock
  _bs->write(_frameWidth/16, 8);
  _bs->write(_frameHeight/16, 8);
  _bs->write(_qp, 8);
  _bs->write(_numFrames, 16);
  _bs->write(_gop, 8);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeWzFrame()
{
  FILE* fReadPtr    = _files->getFile("src")->getFileHandle();

  clock_t timeStart, timeEnd;
  double cpuTime;

  imgpel* currFrame     = _fb->getCurrFrame();
  imgpel* prevFrame     = _fb->getPrevFrame();
  imgpel* currChroma    = _fb->getCurrChroma();
  imgpel* prevChroma    = _fb->getPrevChroma();
  int*    dctFrame      = _fb->getDctFrame();
  int*    quantDctFrame = _fb->getQuantDctFrame();

  int*    residue       = new int[_frameSize];
  int*    chromaResidue = new int[_frameSize>>1];
  int*    dctUFrame     = new int[_frameSize>>2];
  int*    dctVFrame     = new int[_frameSize>>2];
  int*    quantUFrame   = new int[_frameSize>>2];
  int*    quantVFrame   = new int[_frameSize>>2];
  int cw = _frameWidth>>1;
  int ch = _frameHeight>>1;
  int chsize = _frameSize>>2;

  timeStart = clock();

  encodeWzHeader();

  // ---------------------------------------------------------------------
  // Calculate quantization step size
  // ---------------------------------------------------------------------
  computeQuantStep();

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < _numFrames/_gop; keyFrameNo++) {
    fseek(fReadPtr, (3*keyFrameNo*_frameSize)>>1, SEEK_SET);
    fread(prevFrame, _frameSize, 1, fReadPtr);
    fread(prevChroma, _frameSize>>1, 1, fReadPtr);
    for (int idx = 1; idx < _gop; idx++) {
      // Start encoding the WZ frame
      int wzFrameNo = keyFrameNo*_gop + idx;

      cout << "Encoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;

      // Read current frame from the source file
      fseek(fReadPtr, (3*wzFrameNo*_frameSize)>>1, SEEK_SET);
      fread(currFrame, _frameSize, 1, fReadPtr);
      fread(currChroma, _frameSize>>1, 1, fReadPtr);

      // ---------------------------------------------------------------------
      // STAGE 1 - Residual coding & DCT
      // ---------------------------------------------------------------------
      for (int idx = 0; idx < _frameSize; idx++)
        residue[idx] = currFrame[idx] - prevFrame[idx];

      _trans->dctTransform(residue, dctFrame, _frameWidth, _frameHeight);

      updateMaxValue(dctFrame);

      // ---------------------------------------------------------------------
      // STAGE 2 - Quantization
      // ---------------------------------------------------------------------
      _trans->quantization(dctFrame, quantDctFrame, _frameWidth, _frameHeight);

      // ---------------------------------------------------------------------
      // STAGE 3 - Mode decision
      // ---------------------------------------------------------------------
# if MODE_DECISION
      selectCodingMode(quantDctFrame);
# endif // MODE_DECISION

      // ---------------------------------------------------------------------
      // STAGE 4 - Skip mode
      // ---------------------------------------------------------------------
# if SKIP_MODE
      generateSkipMask();

      encodeSkipMask();
# endif // SKIP_MODE

      // ---------------------------------------------------------------------
      // STAGE 5 - Encode (Channel/Entropy)
      // ---------------------------------------------------------------------
      int numBands = 0;

      _rcBitPlaneNum = 0;

      for (int bandNo = 0; bandNo < 16; bandNo++) {
        int x = ScanOrder[bandNo][0];
        int y = ScanOrder[bandNo][1];

        if (bandNo < _numChnCodeBands) {
          _rcQuantMatrix[y][x] = QuantMatrix[_qp][y][x];
          _rcBitPlaneNum += _rcQuantMatrix[y][x];
        }

        if (QuantMatrix[_qp][y][x] > 0)
          numBands++;
      }

      int bits = 0;

      // Entropy encode
      if (numBands > _numChnCodeBands)
        bits = _cavlc->encode(quantDctFrame, _skipMask);

# if HARDWARE_FLOW
      if (bits%32 != 0) {
        int dummy = 32 - (bits%32);
        _bs->write(0, dummy);
      }
# endif // HARDWARE_FLOW

      // Channel encode
      encodeFrameLdpca(quantDctFrame);

      // ---------------------------------------------------------------------
      // STAGE 6 - Write parity and CRC bits to the bitstream
      // ---------------------------------------------------------------------
      for (int i = 0; i < _rcBitPlaneNum; i++)
      {
        for (int j = 0; j < _bitPlaneLength; j++)
          _bs->write(int(_parity[j+i*_bitPlaneLength]), 1);

# if !HARDWARE_FLOW
#   if HARDWARE_LDPC
        if (_bitPlaneLength == 6336)
          for (int n = 0; n < 4; n++)
            _bs->write(_crcPtr[i*4+n], 8);
        else
          _bs->write(_crcPtr[i], 8);
#   else // if !HARDWARE_LDPC
        _bs->write(_crcPtr[i], 8);
#   endif // HARDWARE_LDPC
# endif // !HARDWARE_FLOW
      }
      // Reset crcPtr
      _crcPtr = _crc;

      // ---------------------------------------------------------------------
      // STAGE 7 - Encode Chroma planes and Write to Bit-stream
      // ---------------------------------------------------------------------
      // TODO: implement Chroma encoding
      // - residual enc                       CHECK
      // - DCT (can reuse transforms)         CHECK
      // - Quant                              CHECK  
      // - Entropy Enc                        CHECK
      for (int idx = 0; idx < _frameSize>>1; idx++)
        chromaResidue[idx] = currChroma[idx] - prevChroma[idx];

      memset(dctUFrame, 0, _frameSize>>2);
      memset(dctVFrame, 0, _frameSize>>2);
      memset(quantUFrame, 0, _frameSize>>2);
      memset(quantVFrame, 0, _frameSize>>2);
      _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
      _trans->dctTransform(chromaResidue+chsize, dctVFrame, cw, ch);
      _trans->quantization(dctUFrame, quantUFrame, cw, ch);
      _trans->quantization(dctVFrame, quantVFrame, cw, ch);

      // all bands are CAVLC coded, none channel coded
      _numChnCodeBands = 0;
      int bitsU = _cavlcU->encode(quantUFrame, _bsU);
      int bitsV = _cavlcV->encode(quantVFrame, _bsV);
      cout << (bitsU + bitsV) / 1024  << " Chroma kbits" << endl;
# if HARDWARE_FLOW
      if (bitsU%32 != 0) {
        int dummy = 32 - (bitsU%32);
        _bsU->write(0, dummy);
      }
      if (bitsV%32 != 0) {
        int dummy = 32 - (bitsV%32);
        _bsV->write(0, dummy);
      }
# endif // HARDWARE_FLOW
    } // Finish encoding the WZ frame
  }

  _bs->flush();
  _bsU->flush();
  _bsV->flush();

  timeEnd = clock();
  cpuTime = (timeEnd - timeStart) / CLOCKS_PER_SEC;

  cout << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "Encode statistics" << endl;
  cout << "--------------------------------------------------" << endl;

  report();

  cout << "Total   encoding time: " << cpuTime << "(s)" << endl;
  cout << "Average encoding time: " << cpuTime/_numFrames << "(s)" << endl;
  cout << "--------------------------------------------------" << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::updateMaxValue(int* block)
{
  for (int y = 0; y < 4; y++)
    for (int x = 0; x < 4; x++)
      _maxValue[y][x] = 0;

  for (int y = 0; y < _frameHeight; y += 4)
    for (int x = 0; x < _frameWidth; x += 4)
      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++)
          if (abs(_maxValue[j][i]) < abs(block[(x+i) + (y+j)*_frameWidth]))
            _maxValue[j][i] = block[(x+i) + (y+j)*_frameWidth];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::computeQuantStep()
{
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      if (QuantMatrix[_qp][j][i] != 0) {
        _quantStep[j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
      }
      else
        _quantStep[j][i] = 1;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::selectCodingMode(int* frame)
{
  int numBands = 0;
  int mode;
  int codingMode;

  memset(_average, 0x00, 16*sizeof(double));

  for (int bandNo = 0; bandNo < 16; bandNo++) {
    int i = ScanOrder[bandNo][0];
    int j = ScanOrder[bandNo][1];

    if (QuantMatrix[_qp][j][i] > 0)
      numBands++;
  }

  for (int j = 0; j < _frameHeight; j++)
    for (int i = 0; i < _frameWidth; i++) {
      int mask, data;

      mask = (0x1 << (QuantMatrix[_qp][j%4][i%4]-1)) - 1;
      data = frame[i + j*_frameWidth] & mask;

      _average[(i%4) + (j%4)*4] += data;
    }

  for (int i = 0; i < 16; i++)
    _average[i] /= _bitPlaneLength;

  double th1      = 0.15;
  double th2      = 0.05;
  double th3      = 0.01;
  double energy1  = 0.0;
  double energy2  = 0.0;
  double energy3  = 0.0;

  for (int i = 0; i < 16; i++) {
    int x = ScanOrder[i][0];
    int y = ScanOrder[i][1];

    if (i < 3)
      energy1 += _average[x + y*4];
    else if (i < 6)
      energy2 += _average[x + y*4];
    else
      energy3 += _average[x + y*4];
  }

  energy1 /= 3;

  if (numBands > 3) {
    if (numBands >= 6)
      energy2 /= 3;
    else
      energy2 /= (numBands-3);
  }
  else
    energy2 = 0;

  if (numBands > 6)
    energy3 /= (numBands-6);
  else
    energy3 = 0;

  if (energy1 > (th1/(double)(Scale[0][_qp]))) {
    if (energy2 > (th2/(double)(Scale[1][_qp]))) {
      if (energy3 > (th3/(double)(Scale[2][_qp])))
        mode = 0; // channel coding (channel coding for all bands)
      else
        mode = 2; // hybrid mode 2 (channel coding for lower 6 bands
                        //                entropy coding for other bands)
    }
    else
      mode = 1;   // hybrid mode 1 (channel coding for lower 3 bands
                        //                entropy coding for other bands)
  }
  else
    mode = 3;     // entropy coding (entropy coding for all bands)

  _modeCounter[mode]++;

# if HARDWARE_CMS
  codingMode = _prevMode;
  _prevMode  = mode;
# else // if !HARDWARE_CMS
  codingMode = mode;
# endif // HARDWARE_CMS

  _bs->write(codingMode, 2);

  if (codingMode == 0) _numChnCodeBands = 16; else
  if (codingMode == 1) _numChnCodeBands =  3; else
  if (codingMode == 2) _numChnCodeBands =  6; else
                       _numChnCodeBands =  0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::generateSkipMask()
{
# if HARDWARE_OPT
  int threshold = 2;
# else // if !HARDWARE_OPT
  int threshold = 5;
# endif // HARDWARE_OPT

  int* frame    = _fb->getQuantDctFrame();

  memset(_skipMask, 0x0, sizeof(int)*_bitPlaneLength);

  for (int j = 0; j < _frameHeight; j += SkipBlockSize)
    for (int i = 0; i < _frameWidth; i += SkipBlockSize) {
      int distortion = 0;
      int blockIndex = i/SkipBlockSize + (j/SkipBlockSize)*(_frameWidth/SkipBlockSize);

      for (int y = 0; y < SkipBlockSize; y++)
        for (int x = 0; x < SkipBlockSize; x++) {
          int mask, data;

          mask = (0x1 << (QuantMatrix[_qp][y][x]-1)) - 1;
          data = frame[(i+x) + (j+y)*_frameWidth] & mask;

# if HARDWARE_OPT
          distortion += data;
# else // if !HARDWARE_OPT
          distortion += data * data;
# endif // HARDWARE_OPT
        }

      _skipMask[blockIndex] = (distortion < threshold) ? 1 : 0;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Encoder::encodeSkipMask()
{
  int n0 = 0; // number of non-skipped blocks
  int diff;
  int bitCount;
  int sign;
  int code;
  int length;
  int run = 0;
  int index = 0;

# if HARDWARE_FLOW
  // Pad zero
  int dummy = 20;
  _bs->write(0, dummy);
# endif // HARDWARE_FLOW

  for (int i = 0; i < _bitPlaneLength; i++)
    if (_skipMask[i] == 0)
      n0++;

  diff = abs(n0 - _bitPlaneLength/2);

# if HARDWARE_OPT
  int type = _prevType;
  int nextType = diff / (_bitPlaneLength/6);
  _prevType = nextType;
# else // if !HARDWARE_OPT
  int type = diff / (_bitPlaneLength/6);
# endif // HARDWARE_OPT

  if (type > 2) type = 2;
  _bs->write(type, 2);
  bitCount = 2;

  // Directly output first bit to bitstream
  sign = _skipMask[0];
  _bs->write(sign, 1);
  bitCount++;
  run++;
  index++;

  // Huffman code for other bits
  while (index < _bitPlaneLength) {
    if (_skipMask[index] == sign) {
      run++;
      if (run == 16) { // reach maximum run length
        bitCount += getHuffmanCode(_qp, type, run-1, code, length);
        _bs->write(code, length);
        run = 1;
      }
    }
    else {
      bitCount += getHuffmanCode(_qp, type, run-1, code, length);
      _bs->write(code, length);
      sign = _skipMask[index];
      run = 1;
    }
    index++;
  }

  if (run != 0) {
    bitCount += getHuffmanCode(_qp, type, run-1, code, length);
    _bs->write(code, length);
  }

# if HARDWARE_FLOW
  if (bitCount%32 != 0) { // pad zero
    int dummy = 32 - (bitCount%32);
    _bs->write(0, dummy);
  }
# endif // HARDWARE_FLOW

  return bitCount;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Encoder::getHuffmanCode(int qp, int type, int symbol, int& code, int& length)
{
  int table = qp / 2;

  code   = HuffmanCodeValue [table][type][symbol];
  length = HuffmanCodeLength[table][type][symbol];

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeFrameLdpca(int* frame)
{
  int   bitPosition;
  int*  ldpcaSource = new int[_bitPlaneLength + 8];
  bool* accumulatedSyndrome = _parity;

  for (int band = 0; band < _numChnCodeBands; band++) {
    int i = ScanOrder[band][0];
    int j = ScanOrder[band][1];

    for (bitPosition = _rcQuantMatrix[j][i]-1; bitPosition >= 0; bitPosition--)
    {
      setupLdpcaSource(frame, ldpcaSource, i, j, bitPosition);

# if HARDWARE_LDPC
      if (_bitPlaneLength == 6336) {
        for (int n = 0; n < 4; n++) {
          _ldpca->encode(ldpcaSource + n*1584, accumulatedSyndrome);

          computeCRC(ldpcaSource + n*1584, 1584, _crcPtr+n);

          accumulatedSyndrome += _bitPlaneLength/4;
        }

        _crcPtr += 4;
        cout << ".";
      }
      else {
        _ldpca->encode(ldpcaSource, accumulatedSyndrome);

        cout << ".";

        computeCRC(ldpcaSource, _bitPlaneLength, _crcPtr);

        accumulatedSyndrome += _frameSize/16;

        _crcPtr++;
      }
# else // if !HARDWARE_LDPC
      _ldpca->encode(ldpcaSource, accumulatedSyndrome);

      cout << ".";

      computeCRC(ldpcaSource, _bitPlaneLength, _crcPtr);

      accumulatedSyndrome += _frameSize/16;

      _crcPtr++;
# endif // HARDWARE_LDPC
    }
  }

  cout << endl;

  delete [] ldpcaSource;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::setupLdpcaSource(int* frame, int* source, int offsetX, int offsetY, int bitPosition)
{
  for (int y = 0; y < _frameHeight; y = y+4)
    for (int x = 0; x < _frameWidth; x = x+4) {
      int blockIdx = (x/4) + (y/4)*(_frameWidth/4);
      int frameIdx = (x+offsetX) + (y+offsetY)*_frameWidth;

# if SKIP_MODE

      if (_skipMask[blockIdx] == 1)
        source[blockIdx] = 0;
      else
        source[blockIdx] = (frame[frameIdx] >> bitPosition) & 0x1;

# else // if !SKIP_MODE
      source[blockIdx] = (frame[frameIdx] >> bitPosition) & 0x1;
# endif // SKIP_MODE

    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::computeCRC(int* data, const int length, unsigned char* crc)
{
  // CRC8 110011011
  const int code[9] = {1, 1, 0, 0, 1, 1, 0, 1, 1};

  int* buffer = new int[length + 8];

  memcpy(buffer, data, length*sizeof(int));

  for (int i = length; i < length+8; i++)
    buffer[i] = 0;

  for (int i = 0; i < length; i++)
    if (buffer[i] == 1)
      for (int j = 0; j < 9; j++)
        buffer[i+j] = code[j] ^ buffer[i+j];

  *crc = 0;

  for (int i = 0; i < 8; i++)
    *crc |= buffer[length+i] << (7-i);

  delete [] buffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::report()
{
# if MODE_DECISION
  cout << "Mode usage: ";

  for (int i = 0; i < 4; i++) {
    float usage = (float)_modeCounter[i]/75.0 * 100.0;
    cout << usage << " ";
  }

  cout << endl;
# endif // MODE_DECISION
}
