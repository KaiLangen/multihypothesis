
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "decoder.h"
#include "fileManager.h"
#include "sideInformation.h"
#include "transform.h"
#include "corrModel.h"
#include "time.h"
#include "cavlcDec.h"
#include "frameBuffer.h"
#include "bitstream.h"
#include "ldpcaDec.h"
#include "regExp.h"

using namespace std;


Decoder::Decoder(map<string, string> configMap)
{
  // Parse files
  _files = FileManager::getManager();
  string wzFileName = configMap["WZFile"];
  string recFileName = wzFileName.substr(0, wzFileName.find(".bin")) + ".yuv";
  _files->addFile("wz",     configMap["WZFile"])->openFile("rb");
  _files->addFile("key",    configMap["KeyFile"])->openFile("rb");
  _files->addFile("origin", configMap["SrcFile"])->openFile("rb");
  _files->addFile("rec",    recFileName.c_str())->openFile("wb");

  // create file handles and bitstream Objects for Chroma
  string ubs = wzFileName.substr(0, wzFileName.find(".bin")) + ".u.bin";
  string vbs = wzFileName.substr(0, wzFileName.find(".bin")) + ".v.bin";
  _files->addFile("wzU",  ubs.c_str())->openFile("rb");
  _files->addFile("wzV",  vbs.c_str())->openFile("rb");
  _bs = new Bitstream(1024, _files->getFile("wz")->getFileHandle());
  _bsU = new Bitstream(1024, _files->getFile("wzU")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV")->getFileHandle());

  // Parse other configuration parameters
  _searchParam = atoi(configMap["searchWindowSize"].c_str());
  _searchBlock = atoi(configMap["blockSize"].c_str());

  // compute quantStep
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      if (QuantMatrix[_qp][j][i] != 0) {
        _quantStep[j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
      }
      else
        _quantStep[j][i] = 1;
    }
  }

  decodeWzHeader();

  initialize();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::initialize()
{
  _frameSize        = _frameWidth * _frameHeight;
  _bitPlaneLength   = _frameSize / 16;

  _numChnCodeBands  = 16;

  _dParity          = new double[_bitPlaneLength * BitPlaneNum[_qp]];

# if HARDWARE_LDPC
  if (_bitPlaneLength == 6336)
    _crc            = new unsigned char[BitPlaneNum[_qp] * 4];
  else
    _crc            = new unsigned char[BitPlaneNum[_qp]];
# else
  _crc              = new unsigned char[BitPlaneNum[_qp]];
# endif

  _average          = new double[16];
  _alpha            = new double[_frameSize];
  _sigma            = new double[16];

  _skipMask         = new int[_bitPlaneLength];

  _fb = new FrameBuffer(_frameWidth, _frameHeight, _gop);

  _trans = new Transform(this);

  _model = new CorrModel(this, _trans);
  _si    = new SideInformation(this, _model);

  _cavlc = new CavlcDec(this, 4);

  motionSearchInit(64);

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

  _ldpca = new LdpcaDec(ladderFile, this);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWzHeader()
{
  _frameWidth   = _bs->read(8) * 16;
  _frameHeight  = _bs->read(8) * 16;
  _qp           = _bs->read(8);
  _numFrames    = _bs->read(16);
  _gop          = _bs->read(8);

  cout << "--------------------------------------------------" << endl;
  cout << "WZ frame parameters" << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "Width:  " << _frameWidth << endl;
  cout << "Height: " << _frameHeight << endl;
  cout << "Frames: " << _numFrames << endl;
  cout << "QP:     " << _qp << endl;
  cout << "GOP:    " << _gop << endl;
  cout << "--------------------------------------------------" << endl << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWZframe()
{
  double dPSNRAvg=0;
  double dPSNRSIAvg=0;

  clock_t timeStart, timeEnd;
  double cpuTime;

//  imgpel* currLuma      = _fb->getCurrFrame();
  imgpel* prevLuma      = _fb->getPrevFrame();
  imgpel* currChroma    = _fb->getCurrChroma();
  imgpel* prevChroma    = _fb->getPrevChroma();
  imgpel* oriCurrFrame  = _fb->getorigFrame();
//  imgpel* imgSI         = _fb->getSideInfoFrame();

// Luma Buffers
//  int* iDCT             = _fb->getDctFrame();
//  int* iDCTQ            = _fb->getQuantDctFrame();
//  int* iDecoded         = _fb->getDecFrame();
//  int* iDecodedInvQ     = _fb->getInvQuantDecFrame();
//  int* iDCTBuffer       = new int [_frameSize];
//  int* iDCTResidual     = new int [_frameSize];

// Chroma Buffers
  imgpel* prevKeyLuma   = new imgpel[_frameSize];
  imgpel* prevKeyChroma = new imgpel[_frameSize>>1];
  int* iDecodedU        = new int[_frameSize>>2];
  int* iDecodedV        = new int[_frameSize>>2];
  int* iDctU            = new int[_frameSize>>2];
  int* iDctV            = new int[_frameSize>>2];
  int* iQuantU          = new int[_frameSize>>2];
  int* iQuantV          = new int[_frameSize>>2];

//  int x,y;
  double totalrate=0;
  double dKeyCodingRate=0;
  double dKeyPSNR=0;

  FILE* fReadPtr        = _files->getFile("origin")->getFileHandle();
  FILE* fWritePtr       = _files->getFile("rec")->getFileHandle();
  FILE* fKeyReadPtr     = _files->getFile("key")->getFileHandle();

  parseKeyStat("stats.dat", dKeyCodingRate, dKeyPSNR, _keyQp);

  timeStart = clock();
  int cw = _frameWidth >> 1;
  int ch = _frameWidth >> 1;

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < _numFrames/_gop; keyFrameNo++) {
    // Read previous key frame
    fseek(fKeyReadPtr, (3*(keyFrameNo)*_frameSize)>>1, SEEK_SET);
    fread(prevKeyLuma, _frameSize, 1, fKeyReadPtr);
    fread(prevKeyChroma, _frameSize>>1, 1, fKeyReadPtr);
    fwrite(prevKeyLuma, _frameSize, 1, fWritePtr);
    fwrite(prevKeyChroma, _frameSize>>1, 1, fWritePtr);
    memcpy(prevLuma, prevKeyLuma, _frameSize);
    memcpy(prevChroma, prevKeyChroma, _frameSize>>1);

    for (int idx = 1; idx < _gop; idx++) {
      // Start decoding the WZ frame
      int wzFrameNo = keyFrameNo*_gop + idx;

      cout << "Decoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;
      // Read current frame from the original file (for comparison)
      fseek(fReadPtr, (3*wzFrameNo*_frameSize)>>1, SEEK_SET);
      fread(oriCurrFrame, _frameSize, 1, fReadPtr);

      // ---------------------------------------------------------------------
      // STAGE 1 - Decode Chroma Data
      // ---------------------------------------------------------------------
      int bitsU, bitsV;
      // size of integer buffer in bytes is: (_frameSize >> 2) * 4
      memset(iDecodedU, 0, _frameSize);
      memset(iDecodedV, 0, _frameSize);
      memset(iQuantU, 0, _frameSize);
      memset(iQuantV, 0, _frameSize);

      // read bits from bitstream
      for (int j = 0; j < ch; j += 4)
        for (int i = 0; i < cw; i += 4) {
          bitsU = _cavlc->decode(iDecodedU, i, j, _bsU);
          bitsV = _cavlc->decode(iDecodedV, i, j, _bsV);
        }

# if HARDWARE_FLOW
      if (bitsU%32 != 0) {
        int dummy = 32 - (bitsU%32);
        _bsU->read(dummy);
      }
      if (bitsV%32 != 0) {
        int dummy = 32 - (bitsV%32);
        _bsV->read(dummy);
      }
# endif // HARDWARE_FLOW

      // TODO Fix inverse quantization to prevent seg fault.
      // Maybe create new invQuant function without side-information?
      _trans->invQuantization(iDecodedU, iQuantU, cw, ch);
      _trans->invQuantization(iDecodedV, iQuantV, cw, ch);
      _trans->invDctTransform(iQuantU, iDctU, cw, ch);
      _trans->invDctTransform(iQuantV, iDctV, cw, ch);

      // undo residual encoding
      for (int idx = 0; idx < _frameSize>>2; idx++) {
        currChroma[idx] = iDctU[idx] + prevKeyChroma[idx];
        currChroma[idx+(_frameSize>>2)] = 
          iDctV[idx] + prevKeyChroma[idx+(_frameSize>>2)];
      }
      

      fwrite(oriCurrFrame, _frameSize, 1, fWritePtr);
      fwrite(currChroma, _frameSize>>1, 1, fWritePtr);
      continue;
/*
      // ---------------------------------------------------------------------
      // STAGE 2 - Create side information
      // ---------------------------------------------------------------------
      // Predict from coincident Chroma
      _si->createSideInfo(prevChroma, currChroma, prevLuma, imgSI);

      float currPSNRSI0 = calcPSNR(oriCurrFrame, imgSI, _frameSize);
      cout << "side information quality " << currPSNRSI0 << endl;

      fwrite(imgSI, _frameSize, 1, fWritePtr);
      fwrite(currChroma, _frameSize>>1, 1, fWritePtr);

      // copy curr buffers into prev buffer
      memcpy(prevLuma, imgSI, _frameSize);
      memcpy(prevChroma, currChroma, _frameSize>>1);
      continue;
      // ---------------------------------------------------------------------
      // STAGE 3 - WZ Decode
      // ---------------------------------------------------------------------
      int tmp = getSyndromeData();

      //cout << _numChnCodeBands << endl;

      double dTotalRate = (double)tmp/1024/8;

      _trans->dctTransform(imgSI, iDCT);

      memset(iDecoded, 0, _frameSize*4);
      memset(iDecodedInvQ, 0, _frameSize*4);

      _si->getResidualFrame(prevKeyLuma, imgSI, iDCTBuffer);

      _trans->dctTransform(iDCTBuffer, iDCTResidual);
      _trans->quantization(iDCTResidual, iDCTQ);

      int iOffset = 0;

      for (int i = 0; i < 16; i++) {
        x = ScanOrder[i][0];
        y = ScanOrder[i][1];

#   if MODE_DECISION
        if (i < _numChnCodeBands)
          dTotalRate += decodeLDPC(iDCTQ, iDCTResidual, iDecoded, x, y, iOffset);
#   else
        dTotalRate += decodeLDPC(iDCTQ, iDCTResidual, iDecoded, x, y, iOffset);
#   endif
        iOffset += QuantMatrix[_qp][y][x];
      }

      _trans->invQuantization(iDecoded, iDecodedInvQ, iDCTResidual);
      _trans->invDctTransform(iDecodedInvQ, iDCTBuffer);

      _si->getRecFrame(prevKeyLuma, iDCTBuffer, currLuma);
#     if SKIP_MODE
      getSkippedRecFrame(prevKeyLuma, currLuma, _skipMask);
#     endif

      totalrate += dTotalRate;
      cout << endl;
      cout << "total bytes (Y/frame): " << dTotalRate << " Kbytes" << endl;

      float currPSNRSI = calcPSNR(oriCurrFrame, imgSI, _frameSize);
      cout << "side information quality " << currPSNRSI << endl;
      dPSNRSIAvg += currPSNRSI;

      float currPSNR = calcPSNR(oriCurrFrame, currLuma, _frameSize);
      dPSNRAvg += currPSNR;
      cout << "wyner-ziv frame quality " << currPSNR << endl;
      fwrite(_fb->getCurrFrame(), _frameSize, 1, fWritePtr);
      fwrite(currChroma, _frameSize>>1, 1, fWritePtr);

      // copy curr buffers into prev buffer
      memcpy(prevLuma, currLuma, _frameSize);
      memcpy(prevChroma, currChroma, _frameSize>>1);
*/
    }
  }

  timeEnd = clock();
  cpuTime = (timeEnd - timeStart) / CLOCKS_PER_SEC;

  int iNumGOP = _numFrames/_gop;
  int iDecodeWZFrames = _numFrames - iNumGOP;
  int iTotalFrames = iDecodeWZFrames + iNumGOP;

  cout<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Decode statistics"<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Total Frames        :   "<<iTotalFrames<<endl;
  float framerate = 30.0;
  dPSNRAvg   /= iDecodeWZFrames;
  dPSNRSIAvg /= iDecodeWZFrames;
  cout<<"Total Bytes         :   "<<totalrate<<endl;
  cout<<"WZ Avg Rate  (kbps) :   "<<totalrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0<<endl;
  cout<<"Key Avg Rate (kbps) :   "<<dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Avg Rate (Key+WZ)   :   "<<totalrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0+dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Key Frame Quality   :   "<<dKeyPSNR<<endl;
  cout<<"SI Avg PSNR         :   "<<dPSNRSIAvg<<endl;
  cout<<"WZ Avg PSNR         :   "<<dPSNRAvg<<endl;
  cout<<"Avg    PSNR         :   "<<(dPSNRAvg+dKeyPSNR)/2<<endl;
  cout<<"Total Decoding Time :   "<<cpuTime<<"(s)"<<endl;
  cout<<"Avg Decoding Time   :   "<<cpuTime/(iDecodeWZFrames)<<endl;
  cout<<"--------------------------------------------------"<<endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::parseKeyStat(const char* filename, double &rate, double &psnr, int &QP)
{
  ifstream stats(filename, ios::in);

  if (!stats.is_open())
    return;

  char buf[1024];
  double iSlice_chroma = 0.0;
  double iSliceRate = 0.0;

  RegExp* rgx = RegExp::getInst();

  while (stats.getline(buf, 1024)) {
    string result;

    if (rgx->match(buf, "\\s*SNR Y\\(dB\\)[ |]*([0-9\\.]+)", result)) {
      psnr = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*Average quant[ |]*([0-9\\.]+)", result)) {
      QP = atoi(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*QP[ |]*([0-9\\.]+)", result)) {
      QP = atoi(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*Coeffs\\. C[ |]*([0-9\\.]+)", result)) {
      iSlice_chroma = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*average bits/frame[ |]*([0-9\\.]+)", result)) {
      iSliceRate = atof(result.c_str());
      continue;
    }
  }

  int count = 0;
  double totalRate = 0;

  if (iSliceRate != 0) {
    totalRate += iSliceRate - iSlice_chroma;
    count++;
  }

  if (count)
    rate = totalRate/(1024*count);
  else
    rate = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::getSyndromeData()
{
  int* iDecoded = _fb->getDecFrame();
  int  decodedBits = 0;

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
# if MODE_DECISION
  int codingMode = _bs->read(2);

  if (codingMode == 0) _numChnCodeBands = 16; else
  if (codingMode == 1) _numChnCodeBands =  3; else
  if (codingMode == 2) _numChnCodeBands =  6; else
                       _numChnCodeBands =  0;

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

# else // if !MODE_DECISION

  _rcBitPlaneNum = 0;

  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++) {
      _rcQuantMatrix[j][i] = QuantMatrix[_qp][j][i];
      _rcBitPlaneNum += _rcQuantMatrix[j][i];
    }
# endif // MODE_DECISION

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
# if HARDWARE_FLOW
  // Discard padded bits
  int dummy = 20;
  _bs->read(dummy);

  int bitCount = _bs->getBitCount();
# endif

# if SKIP_MODE
  decodedBits += decodeSkipMask();
# endif

# if HARDWARE_FLOW
  bitCount = _bs->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs->read(dummy);
  }

  bitCount = _bs->getBitCount();
# endif

# if MODE_DECISION
  if (numBands > _numChnCodeBands) {
    for (int j = 0; j < _frameHeight; j += 4)
      for (int i = 0; i < _frameWidth; i += 4) {
        if (_skipMask[i/4+(j/4)*(_frameWidth/4)] == 0) //not skip
          decodedBits += _cavlc->decode(iDecoded, i, j, _bs);
        else
          _cavlc->clearNnz(i/4+(j/4)*(_frameWidth/4));
      }
  }
# endif

# if HARDWARE_FLOW
  bitCount = _bs->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs->read(dummy);
  }
# endif

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Read parity and CRC bits from the bitstream
  for (int i = 0; i < _rcBitPlaneNum; i++)
  {
    for (int j = 0; j < _bitPlaneLength; j++)
      _dParity[j+i*_bitPlaneLength] = (double)_bs->read(1);

# if !HARDWARE_FLOW
#   if HARDWARE_LDPC
    if (_bitPlaneLength == 6336)
      for (int n = 0; n < 4; n++)
        _crc[i*4+n] = _bs->read(8);
    else
      _crc[i] = _bs->read(8);
#   else
    _crc[i] = _bs->read(8);
#   endif
# endif
  }

  return decodedBits;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::decodeSkipMask()
{
  int type   = 0;
  int sign   = 0;
  int index  = 0;
  int length = 0;

  type = _bs->read(2);
  sign = _bs->read(1);

  memset(_skipMask, 0, _bitPlaneLength*sizeof(int));

  while (index < _bitPlaneLength) {
    int code = 0;
    int run  = 0;

    for (length = 1; length < 15; length++) {
      code <<= 1;
      code |= _bs->read(1);

      for (int i = 0; i < 16; i++) {
        int table = _qp / 2;

        if (HuffmanCodeValue [table][type][i] == code &&
            HuffmanCodeLength[table][type][i] == length) {
          run = i+1;
          goto DecodeSkipMaskHuffmanCodeDone;
        }
      }
    }

    DecodeSkipMaskHuffmanCodeDone:

    // Reconstruct skip mask
    if (run == 16)
      for (int i = 0; i < 15; i++)
        _skipMask[index++] = sign;
    else {
      for (int i = 0; i < run; i++)
        _skipMask[index++] = sign;

      sign = (sign == 0) ? 1 : 0;
    }
  }

  return length;
}

/*
*Decoding Process of LDPC
*Param
*/
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double Decoder::decodeLDPC(int* iQuantDCT, int* iDCT, int* iDecoded, int x, int y, int iOffset)
{
  int iCurrPos;
  int* iDecodedTmp;
  double* dLLR;
  double* dAccumulatedSyndrome;
  double* dLDPCDecoded;
  double* dSource;
  double dRate,dTotalRate;
  double dErr;
  double dParityRate;
  unsigned char* ucCRCCode;
  int    iNumCode;

# if HARDWARE_LDPC
  ucCRCCode             = _crc + iOffset*4;
# else
  ucCRCCode             = _crc + iOffset;
# endif
  dAccumulatedSyndrome  = _dParity + _bitPlaneLength*iOffset;

  dLLR         = new double[_bitPlaneLength];
  iDecodedTmp  = new int   [_bitPlaneLength];
  dLDPCDecoded = new double[_bitPlaneLength];
  dSource      = new double[_bitPlaneLength];
  dTotalRate   = 0;

  memset(iDecodedTmp, 0, _bitPlaneLength*4);

  dParityRate = 0;

  for (iCurrPos = _rcQuantMatrix[y][x]-1; iCurrPos >= 0; iCurrPos--)
  {
    if (iCurrPos == _rcQuantMatrix[y][x]-1)
      dParityRate = _model->getSoftInput(iQuantDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 1);
    else
      dParityRate = _model->getSoftInput(iDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 2);
    iNumCode = int(dParityRate*66);

    if (iNumCode <= 2)
      iNumCode = 2;
    if (iNumCode >= 66)
      iNumCode = 66;
    iNumCode = 2;

# if HARDWARE_LDPC
    if (_bitPlaneLength == 6336) {
      double dRateTmp = 0;

      for (int n = 0; n < 4; n++) {
        iNumCode = 2;

        _ldpca->decode(dLLR+n*1584, dAccumulatedSyndrome+n*1584, dSource+n*1584, dLDPCDecoded+n*1584, &dRate, &dErr, *(ucCRCCode+n), iNumCode);
        dRateTmp += (dRate/4.0);
        //cout<<dRate<<endl;
        dRate = 0;
      }
      ucCRCCode += 4;
      cout << ".";

      dTotalRate += dRateTmp;
      dRate = 0;
    }
    else {
      _ldpca->decode(dLLR, dAccumulatedSyndrome, dSource, dLDPCDecoded, &dRate, &dErr, *ucCRCCode, iNumCode);
      cout << ".";

      dTotalRate += dRate;
      dRate = 0;
      ucCRCCode++;
    }
# else
    _ldpca->decode(dLLR, dAccumulatedSyndrome, dSource, dLDPCDecoded, &dRate, &dErr, *ucCRCCode, iNumCode);
    cout << ".";

    dTotalRate += dRate;
    dRate = 0;
    ucCRCCode++;
# endif

    for (int iIndex = 0; iIndex < _bitPlaneLength; iIndex++)
      if (dLDPCDecoded[iIndex] == 1)
        iDecodedTmp[iIndex] |= 0x1<<iCurrPos;

    dAccumulatedSyndrome += _bitPlaneLength;

    memset(dLDPCDecoded, 0, _bitPlaneLength*sizeof(double));
  }

  for (int j = 0; j < _frameHeight; j = j+4)
    for (int i = 0; i < _frameWidth; i = i+4) {
      int tmp = i/4 + j/4*(_frameWidth/4);
      iDecoded[(i+x)+(j+y)*_frameWidth] = iDecodedTmp[tmp];
    }

  delete [] dLLR;
  delete [] iDecodedTmp;
  delete [] dLDPCDecoded;
  delete [] dSource;

  return (dTotalRate*_bitPlaneLength/8/1024);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::getSourceBit(int *dct_q,double *source,int q_i,int q_j,int curr_pos){
  int iWidth,iHeight;
  iWidth  = _frameWidth;
  iHeight = _frameHeight;
  for(int y=0;y<iHeight;y=y+4)
    for(int x=0;x<iWidth;x=x+4)
    {
      source[(x/4)+(y/4)*(iWidth/4)]=(dct_q[(x+q_i)+(y+q_j)*iWidth]>>curr_pos)&(0x1);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::motionSearchInit(int maxsearch_range)
{
  _spiralHpelSearchX = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralSearchX     = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralHpelSearchY = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralSearchY     = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];

  int k,i,l;

  _spiralSearchX[0] = _spiralSearchY[0] = 0;
  _spiralHpelSearchX[0] = _spiralHpelSearchY[0] = 0;

  for (k=1, l=1; l <= std::max<int>(1,maxsearch_range); l++) {
    for (i=-l+1; i< l; i++) {
      _spiralSearchX[k] =  i;
      _spiralSearchY[k] = -l;
      _spiralHpelSearchX[k] =  i<<1;
      _spiralHpelSearchY[k++] = -l<<1;
      _spiralSearchX[k] =  i;
      _spiralSearchY[k] =  l;
      _spiralHpelSearchX[k] =  i<<1;
      _spiralHpelSearchY[k++] =  l<<1;
    }
    for (i=-l;   i<=l; i++) {
      _spiralSearchX[k] = -l;
      _spiralSearchY[k] =  i;
      _spiralHpelSearchX[k] = -l<<1;
      _spiralHpelSearchY[k++] = i<<1;
      _spiralSearchX[k] =  l;
      _spiralSearchY[k] =  i;
      _spiralHpelSearchX[k] =  l<<1;
      _spiralHpelSearchY[k++] = i<<1;
    }
  }
}

void Decoder::getSkippedRecFrame(imgpel* imgPrevKey,imgpel* imgWZFrame, int* skipMask)
{
  for(int j=0;j<_frameHeight;j+=4)
    for(int i=0;i<_frameWidth;i+=4)
    {
      int idx= i/4 + j/4*(_frameWidth/4);
      if(skipMask[idx]==1)//skip
      {
        for(int y=0;y<4;y++)
          for(int x=0;x<4;x++)
          {
            imgWZFrame[(i+x)+(j+y)*_frameWidth]=imgPrevKey[(i+x)+(j+y)*_frameWidth];
          }
      }
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double calcPSNR(unsigned char* img1,unsigned char* img2,int length)
{
  float PSNR;
  float MSE=0;

  for(int i=0;i<length;i++)
    {
      MSE+=pow(float(img1[i]-img2[i]),float(2.0))/length;
    }
  PSNR=10*log10(255*255/MSE);
  //cout<<"PSNR: "<<PSNR<<" dB"<<endl;
  return PSNR;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int getSymbol(int len, int& curr_pos, char* buffer)
{
  int temp = 0;

  for (int count = 0; count < len; count++) {
    int pos = count + curr_pos;

    temp <<= 1;
    temp |= 0x1 & (buffer[pos/8]>>(7-(pos%8)));
  }

  curr_pos += len;

  return temp;
}

