#include <cmath>
#include <cstring>
#include <cassert>
#include <iostream>
#include <map>

#include "calculations.h"
#include "encoder.h"
#include "transform.h"
#include "fileManager.h"
#include "frameBuffer.h"
#include "cavlcEnc.h"
#include "cavlcDec.h"
#include "ldpcaEnc.h"
#include "bitstream.h"

using namespace std;

class EncTest : public Encoder {
public:
  EncTest(map<string, string> configMap);
  ~EncTest() {};
  void resiTest();
  void transTest();
  void qTest1();
  void qTest2();
  void cavlcEnc();
  void cavlcTest1();
  void cavlcTest2();
  void cavlcTestLoop();
  void multiQPTest();

private:
  int* chromaResidue = new int[_frameSize>>1];
  int* iDecodedU = new int[_frameSize>>2];
  int* iDecodedV = new int[_frameSize>>2];
  int* iDctU = new int[_frameSize>>2];
  int* iDctV = new int[_frameSize>>2];
  int* iQuantU = new int[_frameSize>>2];
  int* iQuantV = new int[_frameSize>>2];
  int* rcU = new int[_frameSize/256];
  int* rcV = new int[_frameSize/256];
  int* dctUFrame, *dctVFrame, *quantUFrame, *quantVFrame;
  imgpel* recon = new imgpel[_frameSize>>1];
  const int cw = _frameWidth >> 1;
  const int ch = _frameHeight >> 1;
  const int chsize = _frameSize>>2;
  int bitsU_in;
  int bitsV_in;

  imgpel* currChroma;
  imgpel* prevChroma;
  imgpel* nextChroma;
  FILE* fReadPtr;
};

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EncTest::EncTest(map<string, string> configMap): Encoder(configMap){
  fReadPtr = _files->getFile("src")->getFileHandle();
  currChroma = _fb->getCurrFrame() + _frameSize;
  prevChroma = _fb->getPrevFrame() + _frameSize;
  nextChroma = _fb->getNextFrame() + _frameSize;
  dctUFrame = _fb->getDctFrame() + _frameSize;
  dctVFrame = dctUFrame + chsize;
  quantUFrame = _fb->getQuantDctFrame() + _frameSize;
  quantVFrame = quantUFrame + chsize;

  // for testing file I/O
  _files = FileManager::getManager();
  string ubs = "test.u.bin";
  string vbs = "test.v.bin";
  _files->addFile("wzU_in", ubs.c_str())->openFile("wb");
  _files->addFile("wzV_in", vbs.c_str())->openFile("wb");
  _bsU = new Bitstream(1024, _files->getFile("wzU_in")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV_in")->getFileHandle());
  _qp = 7;
  _chrQp = 7;
  computeQuantStep();

  fseek(fReadPtr, _frameSize, SEEK_SET);
  fread(prevChroma, _frameSize>>1, 1, fReadPtr);
  fseek(fReadPtr, _frameSize, SEEK_CUR);
  fread(currChroma, _frameSize>>1, 1, fReadPtr);
  fseek(fReadPtr, _frameSize, SEEK_CUR);
  fread(nextChroma, _frameSize>>1, 1, fReadPtr);

  rcU = computeResidue(chromaResidue,prevChroma,nextChroma,currChroma,_bsU);
  rcV = computeResidue(chromaResidue+chsize, prevChroma+chsize,
                       nextChroma+chsize, currChroma+chsize, _bsV);
}

// -----------------------------------------------------------------------------
// Test functions
// -----------------------------------------------------------------------------
void EncTest::resiTest() {
  cout << "residual Test...";
  getRecFrame(recon, prevChroma, nextChroma, chromaResidue, rcU, true); 
  getRecFrame(recon+chsize, prevChroma+chsize, nextChroma+chsize,
              chromaResidue+chsize, rcV, true); 

  assert(calcMSE(currChroma, recon, chsize) == 0);
  cout << "passed!" << endl;
}

void EncTest::transTest() {
  cout << "     DCT Test...";
  
  memset(dctUFrame, 0, _frameSize);
  memset(dctVFrame, 0, _frameSize);
  memset(iDctU, 0, _frameSize);
  memset(iDctV, 0, _frameSize);

  // read Chroma from file, then copy into integer buffer
  _trans->dctTransform(chromaResidue, dctUFrame, true);
  _trans->dctTransform(chromaResidue+chsize, dctVFrame, true);
  _trans->invDctTransform(dctUFrame, iDctU, true);
  _trans->invDctTransform(dctVFrame, iDctV, true);

  getRecFrame(recon, prevChroma, nextChroma, iDctU, rcU, true); 
  getRecFrame(recon+chsize, prevChroma+chsize, nextChroma+chsize,iDctV, rcV, true); 

  // compute distortion
  assert(calcMSE(currChroma, recon, chsize) < 1);
  assert(calcMSE(currChroma+chsize, recon+chsize, chsize) < 1);
  cout << "passed!" << endl;
}


void EncTest::qTest1() {
  cout << " Quant Test 1...";
  memset(dctUFrame, 0, _frameSize);
  memset(dctVFrame, 0, _frameSize);
  memset(quantUFrame, 0, _frameSize);
  memset(quantVFrame, 0, _frameSize);
  memset(iDctU, 0, _frameSize);
  memset(iDctV, 0, _frameSize);
  memset(iQuantU, 0, _frameSize);
  memset(iQuantV, 0, _frameSize);

  _trans->dctTransform(chromaResidue, dctUFrame, true);
  _trans->dctTransform(chromaResidue+chsize, dctVFrame, true);
  _trans->quantization(dctUFrame, quantUFrame, true);
  _trans->quantization(dctVFrame, quantVFrame, true);

  _trans->invQuantization(quantUFrame, iQuantU);
  _trans->invQuantization(quantVFrame, iQuantV);

  assert(calcMSE(currChroma, recon, chsize) < 3);
  assert(calcMSE(currChroma+chsize, recon+chsize, chsize) < 3);
  cout << "passed!" << endl;
}

void EncTest::qTest2() {
  cout << " Quant Test 2...";
  memset(dctUFrame, 0, _frameSize);
  memset(dctVFrame, 0, _frameSize);
  memset(quantUFrame, 0, _frameSize);
  memset(quantVFrame, 0, _frameSize);
  memset(iDctU, 0, _frameSize);
  memset(iDctV, 0, _frameSize);
  memset(iQuantU, 0, _frameSize);
  memset(iQuantV, 0, _frameSize);

  _trans->dctTransform(chromaResidue, dctUFrame, true);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, true);
  _trans->quantization(dctUFrame, quantUFrame, true);
  _trans->quantization(dctVFrame, quantVFrame, true);

  _trans->invQuantization(quantUFrame, iQuantU);
  _trans->invQuantization(quantVFrame, iQuantV);
  _trans->invDctTransform(iQuantU, iDctU, true);
  _trans->invDctTransform(iQuantV, iDctV, true);

  getRecFrame(recon, prevChroma, nextChroma, iDctU, rcU, true); 
  getRecFrame(recon+chsize, prevChroma+chsize, nextChroma+chsize,iDctV, rcV, true); 

  assert(calcMSE(currChroma, recon, chsize) < 3);
  assert(calcMSE(currChroma+chsize, recon+chsize, chsize) < 3);
  cout << "passed!" << endl;
}

void EncTest::cavlcEnc() {
  cout << " Cavlc Encode...";

  memset(dctUFrame, 0, _frameSize);
  memset(dctVFrame, 0, _frameSize);
  memset(quantUFrame, 0, _frameSize);
  memset(quantVFrame, 0, _frameSize);

  _trans->dctTransform(chromaResidue, dctUFrame, true);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, true);
  _trans->quantization(dctUFrame, quantUFrame, true);
  _trans->quantization(dctVFrame, quantVFrame, true);
  _numChnCodeBands = 0;
  CavlcEnc* cavlcEncU = new CavlcEnc(this, 4);
  CavlcEnc* cavlcEncV = new CavlcEnc(this, 4);
  bitsU_in = cavlcEncU->encode(quantUFrame, _bsU);
  bitsV_in = cavlcEncV->encode(quantVFrame, _bsV);
  _bsU->flush();
  _bsV->flush();
  _files->getFile("wzU_in")->closeFile();
  _files->getFile("wzV_in")->closeFile();
  cout << "passed!" << endl;
}

void EncTest::cavlcTest1() {
  cout << " Cavlc Test 1...";

  memset(dctUFrame, 0, _frameSize);
  memset(dctVFrame, 0, _frameSize);
  memset(quantUFrame, 0, _frameSize);
  memset(quantVFrame, 0, _frameSize);
  memset(iDecodedU, 0, _frameSize);
  memset(iDecodedV, 0, _frameSize);
  string ubs = "test.u.bin";
  string vbs = "test.v.bin";
  _files->addFile("wzU_out", ubs.c_str())->openFile("rb");
  _files->addFile("wzV_out", vbs.c_str())->openFile("rb");
  Bitstream* bsU_out = new Bitstream(1024, _files->getFile("wzU_out")->getFileHandle());
  Bitstream* bsV_out = new Bitstream(1024, _files->getFile("wzV_out")->getFileHandle());
  CavlcDec* cavlcDecU = new CavlcDec(this, 4);
  CavlcDec* cavlcDecV = new CavlcDec(this, 4);

  // Bits may have been written to stream as part of rcList encoding
  // read the rcList from Bitstreams 
# if !HARDWARE_FLOW
  assert(rcU);
  assert(rcV);
  for (int i = 0; i < _frameSize/256; i++) {
    assert(bsU_out->read(1) == rcU[i]);
    assert(bsV_out->read(1) == rcV[i]);
  }
#endif

  _trans->dctTransform(chromaResidue, dctUFrame, true);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, true);
  _trans->quantization(dctUFrame, quantUFrame, true);
  _trans->quantization(dctVFrame, quantVFrame, true);

  int bitsU_out = 0;
  int bitsV_out = 0;
  // read bits from bitstream
  for (int j = 0; j < ch; j += 4) {
    for (int i = 0; i < cw; i += 4) {
      bitsU_out += cavlcDecU->decode(iDecodedU, i, j, bsU_out);
      bitsV_out += cavlcDecV->decode(iDecodedV, i, j, bsV_out);
    }
  }
  assert(bitsU_in == bitsU_out);
  assert(bitsV_in == bitsV_out);
  assert(calcMSE(quantUFrame, iDecodedU, chsize) == 0);
  assert(calcMSE(quantVFrame, iDecodedV, chsize) == 0);
  cout << "passed!" << endl;
}

void EncTest::cavlcTest2() {
  cout << " Cavlc Test 2...";
  memset(iDecodedU, 0, _frameSize);
  memset(iDecodedV, 0, _frameSize);
  memset(iDctU, 0, _frameSize);
  memset(iDctV, 0, _frameSize);
  memset(iQuantU, 0, _frameSize);
  memset(iQuantV, 0, _frameSize);
  string ubs = "test.u.bin";
  string vbs = "test.v.bin";
  _files->addFile("wzU_out", ubs.c_str())->openFile("rb");
  _files->addFile("wzV_out", vbs.c_str())->openFile("rb");
  Bitstream* bsU_out = new Bitstream(1024, _files->getFile("wzU_out")->getFileHandle());
  Bitstream* bsV_out = new Bitstream(1024, _files->getFile("wzV_out")->getFileHandle());
  CavlcDec* cavlcDecU = new CavlcDec(this, 4);
  CavlcDec* cavlcDecV = new CavlcDec(this, 4);

  // Bits may have been written to stream as part of rcList encoding
  // read the rcList from Bitstreams 
# if !HARDWARE_FLOW
  assert(rcU);
  assert(rcV);
  for (int i = 0; i < _frameSize/256; i++) {
    assert(bsU_out->read(1) == rcU[i]);
    assert(bsV_out->read(1) == rcV[i]);
  }
#endif

  // cavlc decode
  int bitsU_out = 0;
  int bitsV_out = 0;
  for (int j = 0; j < ch; j += 4) {
    for (int i = 0; i < cw; i += 4) {
      bitsU_out += cavlcDecU->decode(iDecodedU, i, j, bsU_out);
      bitsV_out += cavlcDecV->decode(iDecodedV, i, j, bsV_out);
    }
  }

  // inverse Q & DCT
  _trans->invQuantization(iDecodedU, iQuantU);
  _trans->invQuantization(iDecodedV, iQuantV);
  _trans->invDctTransform(iQuantU, iDctU, true);
  _trans->invDctTransform(iQuantV, iDctV, true);

  getRecFrame(recon, prevChroma, nextChroma, iDctU, rcU, true); 
  getRecFrame(recon+chsize, prevChroma+chsize, nextChroma+chsize,iDctV, rcV, true); 

  assert(calcMSE(currChroma, recon, chsize) < 3);
  assert(calcMSE(currChroma+chsize, recon+chsize, chsize) < 3);
  cout << "passed!" << endl;
}


void EncTest::cavlcTestLoop() {
  int minf = 5;
  int maxf = 9;
  cout << "####################################" << endl; 
  cout << " Residual encode/decode loop" << endl;
  cout << "####################################" << endl; 
  bitsU_in = 0;
  bitsV_in = 0;
  _files->getFile("wzU_in")->openFile("wb");
  _files->getFile("wzV_in")->openFile("wb");
  _bsU = new Bitstream(1024, _files->getFile("wzU_in")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV_in")->getFileHandle());
  CavlcEnc* cavlcEncU = new CavlcEnc(this, 4);
  CavlcEnc* cavlcEncV = new CavlcEnc(this, 4);
  _numChnCodeBands = 0;
  cout << " Encode:" << endl;
  for(int frIdx = minf; frIdx < maxf; frIdx++) {
    memset(dctUFrame, 0, _frameSize);
    memset(dctVFrame, 0, _frameSize);
    memset(quantUFrame, 0, _frameSize);
    memset(quantVFrame, 0, _frameSize);

    fseek(fReadPtr, frIdx*(3*_frameSize>>1), SEEK_SET);
    fseek(fReadPtr, _frameSize, SEEK_CUR);
    fread(currChroma, _frameSize>>1, 1, fReadPtr);

    for (int idx = 0; idx < _frameSize>>1; idx++)
      chromaResidue[idx] = currChroma[idx] - prevChroma[idx];

    _trans->dctTransform(chromaResidue, dctUFrame, true);
    _trans->dctTransform(chromaResidue+chsize, dctVFrame, true);
    _trans->quantization(dctUFrame, quantUFrame, true);
    _trans->quantization(dctVFrame, quantVFrame, true);
    int tempU = cavlcEncU->encode(quantUFrame, _bsU);
    int tempV = cavlcEncV->encode(quantVFrame, _bsV);
    cout << "Frame " << frIdx << " Bits (U): " << tempU << endl;
    cout << "Frame " << frIdx << " Bits (V): " << tempV << endl;
    bitsU_in += tempU;
    bitsV_in += tempV;
  }
  cout << "Bits written (U): " << bitsU_in << endl;
  cout << "Bits written (V): " << bitsV_in << endl << endl;
  _bsU->flush();
  _bsV->flush();
  _files->getFile("wzU_in")->closeFile();
  _files->getFile("wzV_in")->closeFile();

  /*********************************************************/

  _files->getFile("wzU_out")->openFile("rb");
  _files->getFile("wzV_out")->openFile("rb");
  Bitstream* bsU_out = new Bitstream(1024, _files->getFile("wzU_out")->getFileHandle());
  Bitstream* bsV_out = new Bitstream(1024, _files->getFile("wzV_out")->getFileHandle());
  CavlcDec* cavlcDecU = new CavlcDec(this, 4);
  CavlcDec* cavlcDecV = new CavlcDec(this, 4);
  int bitsU_out = 0;
  int bitsV_out = 0;
  cout << " Decode:" << endl;
  for(int frIdx = minf; frIdx < maxf; frIdx++) {
    memset(iDecodedU, 0, _frameSize);
    memset(iDecodedV, 0, _frameSize);
    memset(iDctU, 0, _frameSize);
    memset(iDctV, 0, _frameSize);
    memset(iQuantU, 0, _frameSize);
    memset(iQuantV, 0, _frameSize);

    // cavlc decode
    int tempU = 0;
    int tempV = 0;
    for (int j = 0; j < ch; j += 4) {
      for (int i = 0; i < cw; i += 4) {
        tempU += cavlcDecU->decode(iDecodedU, i, j, bsU_out);
        tempV += cavlcDecV->decode(iDecodedV, i, j, bsV_out);
      }
    }
    bitsU_out += tempU;
    bitsV_out += tempV;
    cout << "Frame " << frIdx << " Bits (U): " << tempU << endl;
    cout << "Frame " << frIdx << " Bits (V): " << tempV << endl;

    // inverse Q & DCT
    _trans->invQuantization(iDecodedU, iQuantU);
    _trans->invQuantization(iDecodedV, iQuantV);
    _trans->invDctTransform(iQuantU, iDctU, true);
    _trans->invDctTransform(iQuantV, iDctV, true);

    for (int idx = 0; idx < chsize; idx++) {
      recon[idx] = iDctU[idx] + prevChroma[idx];
      recon[idx+chsize] = iDctV[idx] + prevChroma[idx+chsize];
    }

    double currPSNR = calcPSNR(currChroma, recon, chsize);
    double currMSE = calcMSE(currChroma, recon, chsize);
    cout << "PSNR (U): " << currPSNR << endl;
    cout << "MSE (U): " << currMSE << endl;
    currPSNR = calcPSNR(currChroma+chsize, recon+chsize, chsize);
    currMSE = calcMSE(currChroma+chsize, recon+chsize, chsize);
    cout << "PSNR (V): " << currPSNR << endl;
    cout << "MSE (V): " << currMSE << endl;
  }
}

void EncTest::multiQPTest() {
  cout << "####################################" << endl; 
  cout << " Residual encode/decode loop" << endl;
  cout << "####################################" << endl; 
  fseek(fReadPtr, _frameSize, SEEK_CUR);
  fread(currChroma, _frameSize>>1, 1, fReadPtr);
  for (int idx = 0; idx < _frameSize>>1; idx++)
    chromaResidue[idx] = currChroma[idx] - prevChroma[idx];

  for(int qp = 7; qp >= 0; qp--) {
    _chrQp = qp;
    computeQuantStep();
    memset(dctUFrame, 0, _frameSize);
    memset(dctVFrame, 0, _frameSize);
    memset(quantUFrame, 0, _frameSize);
    memset(quantVFrame, 0, _frameSize);
    memset(iDecodedU, 0, _frameSize);
    memset(iDecodedV, 0, _frameSize);
    memset(iDctU, 0, _frameSize);
    memset(iDctV, 0, _frameSize);
    memset(iQuantU, 0, _frameSize);
    memset(iQuantV, 0, _frameSize);

    cout << " Encode:" << endl;
    _trans->dctTransform(chromaResidue, dctUFrame, true);
    _trans->dctTransform(chromaResidue+chsize, dctVFrame, true);
    _trans->quantization(dctUFrame, quantUFrame, true);
    _trans->quantization(dctVFrame, quantVFrame, true);

    cout << " Decode:" << endl;
    _trans->invQuantization(quantUFrame, iQuantU);
    _trans->invQuantization(quantVFrame, iQuantV);
    double currPSNR = calcPSNR(chromaResidue, iDctU, chsize);
    double currMSE = calcMSE(chromaResidue, iDctU, chsize);
    _trans->invDctTransform(iQuantU, iDctU, true);
    _trans->invDctTransform(iQuantV, iDctV, true);

    currPSNR = calcPSNR(chromaResidue, iDctU, chsize);
    currMSE = calcMSE(chromaResidue, iDctU, chsize);
    cout <<  "MSE (U): " << currMSE << endl;
    cout << "PSNR (U): " << currPSNR << endl;
    currPSNR = calcPSNR(chromaResidue+chsize, iDctV, chsize);
    currMSE = calcMSE(chromaResidue+chsize, iDctV, chsize);
    cout << "MSE (V): " << currMSE << endl;
    cout << "PSNR (V): " << currPSNR << endl << endl;
  }
  _chrQp = 7;
  computeQuantStep();
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc != 3) {
    cerr << endl;
    cerr << "Usage: ./tester ";
    cerr << "[config file] ";
    cerr << "[input video]";
    cerr << endl;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1], true);
    configMap["SrcFile"] = argv[2];
    EncTest* test = new EncTest(configMap);

    test->resiTest();
    test->transTest();
    test->qTest1();
    test->qTest2();
//    test->multiQPTest();

    test->cavlcEnc();
    test->cavlcTest1();
    test->cavlcTest2();
//    test->cavlcTestLoop();

  }
  return 0;
}
