
#include <iostream>
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

#include "decoder.h"
#include "fileManager.h"
#include "sideInformation.h"
#include "transform.h"
#include "corrModel.h"
#include "calculations.h"
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

  // read height, width, gop, and qp's from header
  decodeWzHeader();

  _files->addFile("oracle", configMap["OracleFile"])->openFile("rb");

  // Parse other configuration parameters
  int nRefFrames = atoi(configMap["NumRefFrames"].c_str());
  _trans = new Transform(this);
  _model = new CorrModel(this, _trans);
  _si    = new SideInformation(this, _model, configMap);
  _fb = new FrameBuffer(_frameWidth, _frameHeight, nRefFrames);
  _frameSize        = _frameWidth * _frameHeight;
  _bitPlaneLength   = _frameSize / 16;
  motionSearchInit(64);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWzHeader()
{
  _frameWidth   = _bs->read(8) * 16;
  _frameHeight  = _bs->read(8) * 16;
  _qp           = _bs->read(8);
  _chrQp        = _bs->read(8);
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
void Decoder::decodeWzFrame()
{
  double dPSNRAvg=0;
  double dPSNRUAvg=0;
  double dPSNRVAvg=0;
  double dPSNRSIAvg=0;
  double dPSNRPrevSIAvg=0;
  double dPSNRPrevUAvg=0;
  double dPSNRPrevVAvg=0;
  int iDecodeWZFrames=0;

  clock_t timeStart, timeEnd;
  double cpuTime;
  int fullFrameSize = 3*(_frameSize>>1);
  int chsize = _frameSize>>2;

// Luma Buffers
  imgpel* oriCurrFrame  = _fb->getorigFrame();
  imgpel* currFrame     = _fb->getCurrFrame();
  imgpel* prevKey       = _fb->getPrevFrame();
  imgpel* nextKey       = _fb->getNextFrame();
  imgpel* imgSI         = _fb->getSideInfoFrame();
  imgpel* nextFrame     = new imgpel[fullFrameSize];
  imgpel* prevFrame     = new imgpel[fullFrameSize];

  RefBuffer* refFrames  = _fb->getRefBuffer();

  double totalrate=0;
  double chromarate=0;
  double dKeyCodingRate[3] = {0., 0., 0.};
  double dKeyPSNR[3] = {0., 0., 0.};

  FILE* fReadPtr        = _files->getFile("origin")->getFileHandle();
  FILE* fWritePtr       = _files->getFile("rec")->getFileHandle();
  FILE* fKeyReadPtr     = _files->getFile("key")->getFileHandle();
  FILE* oracleReadPtr = NULL;
  oracleReadPtr = _files->getFile("oracle")->getFileHandle();

  parseKeyStat("stats.dat", dKeyCodingRate, dKeyPSNR);

  timeStart = clock();

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < _numFrames/_gop; keyFrameNo++) {
    // Read previous and next key frame
    fseek(fKeyReadPtr, (3*(keyFrameNo)*_frameSize)>>1, SEEK_SET);
    fread(prevKey, fullFrameSize, 1, fKeyReadPtr);
    fread(nextKey, fullFrameSize, 1, fKeyReadPtr);

    // write prevKey to output file
    fwrite(prevKey, fullFrameSize, 1, fWritePtr);

    // copy into frame buffers
    memcpy(prevFrame, prevKey, fullFrameSize);
    memcpy(nextFrame, nextKey, fullFrameSize);

    refFrames->initPrevNextBuffers();
    int idx = 1;
    while (idx < _gop) {
      // Start decoding the WZ frame
      int wzFrameNo = keyFrameNo*_gop + idx;

      cout << "Decoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;
      // Read current frame from the original file (for comparison)
      fseek(fReadPtr, (3*wzFrameNo*_frameSize)>>1, SEEK_SET);
      fread(oriCurrFrame, fullFrameSize, 1, fReadPtr);

      // ---------------------------------------------------------------------
      // STAGE 2 - Create side information
      // ---------------------------------------------------------------------
      // Predict from coincident Chroma
       fseek(oracleReadPtr, (3*(wzFrameNo)*_frameSize)>>1, SEEK_SET);
       fread(currFrame, fullFrameSize, 1, oracleReadPtr);
       _si->oracle_MEMC(refFrames, imgSI);

      float currPSNRU = calcPSNR(oriCurrFrame+_frameSize, imgSI, chsize);
      cout << "PSNR Recoloured Chroma (U): " << currPSNRU << endl;

      float currPSNRV = calcPSNR(oriCurrFrame+_frameSize+chsize, imgSI+chsize, chsize);
      cout << "PSNR Recoloured Chroma (V): " << currPSNRU << endl;
      dPSNRSIAvg += (currPSNRU + currPSNRV) / 2;
      memcpy(currFrame+_frameSize, imgSI, _frameSize>>1);
//      dPSNRAvg += calcPSNR(oriCurrFrame, currFrame, fullFrameSize);

      float currPSNRLuma = calcPSNR(oriCurrFrame, currFrame, _frameSize);
      cout << "PSNR Luma:" << currPSNRLuma << endl;
      cout << "PSNR Frame Avg: ";
      cout << (currPSNRLuma*6 + currPSNRU + currPSNRV) / 8 << endl << endl;

      // update frameBuffer
      refFrames->updateRecWindow();
      
//      // SI was generated using Chroma-ME, go back a frame
//      if (idx % 2 == 0) {
//        memcpy(nextFrame, currFrame, fullFrameSize);
//        idx--;
//      }
//      // SI was generated using MCI
//      else if (idx < _gop-1) {
//        // update prevFrame ptr
//        memcpy(prevFrame, nextFrame, fullFrameSize);
//
//        // write curr and next frames, then skip forward three frame
//        fwrite(currFrame, fullFrameSize, 1, fWritePtr);
//        fwrite(nextFrame, fullFrameSize, 1, fWritePtr);
//        idx += 3;
//      }
//      // SI was generated using MCI and curr frame is LAST in the GOP
//      else {}
        // write only the current frame to output
        fwrite(currFrame, fullFrameSize, 1, fWritePtr);
        idx++;
        iDecodeWZFrames++;
    }
  }

  timeEnd = clock();
  cpuTime = (timeEnd - timeStart) / CLOCKS_PER_SEC;

  int iNumGOP = _numFrames/_gop;
  int iTotalFrames = iDecodeWZFrames + iNumGOP;

  cout<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Decode statistics"<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Total Frames        :   "<<iTotalFrames<<endl;
  float framerate = 30.0;
  cout<<"Total KBytes        :   "<<totalrate<<endl;
  cout<<"Total Decoding Time :   "<<cpuTime<<"(s)"<<endl;
  cout<<"Avg Decoding Time   :   "<<cpuTime/(iDecodeWZFrames)<<endl;
  cout<<"Y Key Rate (kbps)   :   ";
  cout<<dKeyCodingRate[0]*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"U Key Rate (kbps)   :   ";
  cout<<dKeyCodingRate[2]*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"V Key Rate (kbps)   :   ";
  cout<<dKeyCodingRate[2]*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Key Frame PSNRY     :   "<<dKeyPSNR[0]<<endl;
  cout<<"Key Frame PSNRU     :   "<<dKeyPSNR[1]<<endl;
  cout<<"Key Frame PSNRV     :   "<<dKeyPSNR[2]<<endl;
  dPSNRAvg   /= iDecodeWZFrames;
  cout<<"Avg Rate (Key+WZ)   :   "<<totalrate*framerate/
                                   (double)iTotalFrames*8.0+
                                   dKeyCodingRate[0]*framerate*
                                   (iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Avg    PSNR         :   "<<(dPSNRAvg+dKeyPSNR[0])/2<<endl;
  cout<<"--------------------------------------------------"<<endl;
  dPSNRSIAvg /= iDecodeWZFrames;
  dPSNRUAvg /= iDecodeWZFrames;
  dPSNRVAvg /= iDecodeWZFrames;
  cout<<"PROPOSED SOLUTION   :"<< endl;
  cout<<"Luma SI Avg PSNR    :   "<<dPSNRSIAvg<<endl;
  cout<<"Chroma (U) Avg PSNR :   "<<dPSNRUAvg<<endl;
  cout<<"Chroma (V) Avg PSNR :   "<<dPSNRVAvg<<endl;
  cout<<"Chroma rate (kbps)  :   "<<chromarate*(double)framerate/(double)iTotalFrames<<endl;
  cout<<"--------------------------------------------------"<<endl;
  dPSNRPrevSIAvg /= iDecodeWZFrames;
  dPSNRPrevUAvg /= iDecodeWZFrames;
  dPSNRPrevVAvg /= iDecodeWZFrames;
  cout<<"COPY PREV FRAME     :"<< endl;
  cout<<"Luma SI Avg PSNR    :   "<<dPSNRPrevSIAvg<<endl;
  cout<<"Chroma (U) Avg PSNR :   "<<dPSNRPrevUAvg<<endl;
  cout<<"Chroma (V) Avg PSNR :   "<<dPSNRPrevVAvg<<endl;
  cout<<"--------------------------------------------------"<<endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::parseKeyStat(const char* filename, double* rate, double* psnr)
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
      psnr[0] = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*SNR U\\(dB\\)[ |]*([0-9\\.]+)", result)) {
      psnr[1] = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*SNR V\\(dB\\)[ |]*([0-9\\.]+)", result)) {
      psnr[2] = atof(result.c_str());
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

  double totalRate = 0;

  if (iSliceRate != 0) {
    totalRate += iSliceRate - iSlice_chroma;
    rate[0] = totalRate/(1000);
    rate[1] = iSlice_chroma/(2*1000);
    rate[2] = iSlice_chroma/(2*1000);
  } else {
    rate[0] = 0;
    rate[1] = 0;
    rate[2] = 0;
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
