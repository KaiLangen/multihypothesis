#include <cmath>
#include <cstring>
#include <cassert>
#include <iostream>
#include <map>
#include <opencv2/opencv.hpp>

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
  void transTest();
  void qTest1();
  void qTest2();
  void cavlcEnc();
  void cavlcTest1();
  void cavlcTest2();

private:
  int* chromaResidue = new int[_frameSize>>1];
  int* dctUFrame = new int[_frameSize>>2];
  int* dctVFrame = new int[_frameSize>>2];
  int* quantUFrame = new int[_frameSize>>2];
  int* quantVFrame = new int[_frameSize>>2];
  int* iDecodedU = new int[_frameSize>>2];
  int* iDecodedV = new int[_frameSize>>2];
  int* iDctU = new int[_frameSize>>2];
  int* iDctV = new int[_frameSize>>2];
  int* iQuantU = new int[_frameSize>>2];
  int* iQuantV = new int[_frameSize>>2];
  int* skipMask = new int[_bitPlaneLength];
  imgpel* recon = new imgpel[_frameSize>>1];
  const int cw = _frameWidth >> 1;
  const int ch = _frameHeight >> 1;
  const int chsize = _frameSize>>2;
  int bitsU_in;
  int bitsV_in;

  imgpel* currChroma;
  imgpel* prevChroma;
  FILE* fReadPtr;
};

// -----------------------------------------------------------------------------
//  Helper / Constructor Functions
// -----------------------------------------------------------------------------
template <class T, class U>
double calcMSE(T* img1, U* img2,int length)
{
  float MSE=0;
  for(int i=0;i<length;i++)
  {
    MSE+=pow(float(img1[i]-img2[i]),float(2.0))/length;
  }
  return MSE;
}

template <class T, class U>
double calcPSNR(T* img1, U* img2,int length)
{
  return 10*log10(255*255/calcMSE(img1, img2, length));
}

EncTest::EncTest(map<string, string> configMap): Encoder(configMap){
  fReadPtr = _files->getFile("src")->getFileHandle();
  prevChroma = _fb->getPrevChroma();
  currChroma = _fb->getCurrChroma();

  // for testing file I/O
  _files = FileManager::getManager();
  string ubs = "test.u.bin";
  string vbs = "test.v.bin";
  _files->addFile("wzU_in", ubs.c_str())->openFile("wb");
  _files->addFile("wzV_in", vbs.c_str())->openFile("wb");
  _bsU = new Bitstream(1024, _files->getFile("wzU_in")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV_in")->getFileHandle());
  _qp = 7;
  computeQuantStep();

  fseek(fReadPtr, _frameSize, SEEK_SET);
  fread(prevChroma, _frameSize>>1, 1, fReadPtr);
  fseek(fReadPtr, (15*_frameSize)>>1, SEEK_CUR);
  fseek(fReadPtr, _frameSize, SEEK_CUR);
  fread(currChroma, _frameSize>>1, 1, fReadPtr);
  for (int idx = 0; idx < _frameSize>>1; idx++)
    chromaResidue[idx] = currChroma[idx];
}

// -----------------------------------------------------------------------------
// Test functions
// -----------------------------------------------------------------------------
void EncTest::transTest() {
  cout << "#############################" << endl; 
  cout << "dct->idct:" << endl;
  cout << "#############################" << endl; 
  cout << "Original vs. reconstructed:" << endl;
  
  memset(dctUFrame, 0, _frameSize>>2);
  memset(dctVFrame, 0, _frameSize>>2);
  memset(iDctU, 0, _frameSize>>2);
  memset(iDctV, 0, _frameSize>>2);

  // read Chroma from file, then copy into integer buffer
  _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, cw, ch);
  _trans->invDctTransform(dctUFrame, iDctU, cw, ch);
  _trans->invDctTransform(dctUFrame, iDctU, cw, ch);
  _trans->invDctTransform(dctVFrame, iDctV, cw, ch);

  // compute distortion
  double currPSNR = calcPSNR(currChroma, iDctU, chsize);
  double currMSE = calcMSE(currChroma, iDctU, chsize);
  cout <<  "MSE (U): " << currMSE << endl;
  cout << "PSNR (U): " << currPSNR << endl;
  currPSNR = calcPSNR(currChroma + chsize, iDctV, chsize);
  currMSE = calcMSE(currChroma + chsize, iDctV, chsize);
  cout << "MSE (V): " << currMSE << endl;
  cout << "PSNR (V): " << currPSNR << endl << endl;
}


void EncTest::qTest1() {
  cout << "#############################" << endl; 
  cout << "dct->quant->iquant" << endl;
  cout << "#############################" << endl; 
  cout << "DCT vs. un-quantized:" << endl;
  memset(dctUFrame, 0, _frameSize>>2);
  memset(dctVFrame, 0, _frameSize>>2);
  memset(quantUFrame, 0, _frameSize>>2);
  memset(quantVFrame, 0, _frameSize>>2);
  memset(iDctU, 0, _frameSize>>2);
  memset(iDctV, 0, _frameSize>>2);
  memset(iQuantU, 0, _frameSize>>2);
  memset(iQuantV, 0, _frameSize>>2);

  _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, cw, ch);
  _trans->quantization(dctUFrame, quantUFrame, cw, ch);
  _trans->quantization(dctVFrame, quantVFrame, cw, ch);

  _trans->invQuantization(quantUFrame, iQuantU, cw, ch);
  _trans->invQuantization(quantVFrame, iQuantV, cw, ch);

  // compute distortion
  double currPSNR = calcPSNR(dctUFrame, iQuantU, chsize);
  double currMSE = calcMSE(dctUFrame, iQuantU, chsize);
  cout <<  "MSE (U): " << currMSE << endl;
  cout << "PSNR (U): " << currPSNR << endl;
  currPSNR = calcPSNR(dctVFrame, iQuantV, chsize);
  currMSE = calcMSE(dctVFrame, iQuantV, chsize);
  cout << "MSE (V): " << currMSE << endl;
  cout << "PSNR (V): " << currPSNR << endl << endl;
}


void EncTest::qTest2() {
  cout << "#############################" << endl; 
  cout << "dct->quant->iquant->idct:" << endl;
  cout << "#############################" << endl; 
  cout << "Original vs. reconstructed:" << endl;
  memset(dctUFrame, 0, _frameSize>>2);
  memset(dctVFrame, 0, _frameSize>>2);
  memset(quantUFrame, 0, _frameSize>>2);
  memset(quantVFrame, 0, _frameSize>>2);
  memset(iDctU, 0, _frameSize>>2);
  memset(iDctV, 0, _frameSize>>2);
  memset(iQuantU, 0, _frameSize>>2);
  memset(iQuantV, 0, _frameSize>>2);

  _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, cw, ch);
  _trans->quantization(dctUFrame, quantUFrame, cw, ch);
  _trans->quantization(dctVFrame, quantVFrame, cw, ch);

  _trans->invQuantization(quantUFrame, iQuantU, cw, ch);
  _trans->invQuantization(quantVFrame, iQuantV, cw, ch);
  _trans->invDctTransform(iQuantU, iDctU, cw, ch);
  _trans->invDctTransform(iQuantV, iDctV, cw, ch);

  // compute distortion
  double currPSNR = calcPSNR(currChroma, iDctU, chsize);
  double currMSE = calcMSE(currChroma, iDctU, chsize);
  cout <<  "MSE (U): " << currMSE << endl;
  cout << "PSNR (U): " << currPSNR << endl;
  currPSNR = calcPSNR(currChroma+chsize, iDctV, chsize);
  currMSE = calcMSE(currChroma+chsize, iDctV, chsize);
  cout << "MSE (V): " << currMSE << endl;
  cout << "PSNR (V): " << currPSNR << endl << endl;
}

void EncTest::cavlcEnc() {
  cout << "####################################" << endl; 
  cout << "dct->quant->cavlcEnc->FILE:" << endl;
  cout << "####################################" << endl; 
  memset(dctUFrame, 0, _frameSize>>2);
  memset(dctVFrame, 0, _frameSize>>2);
  memset(quantUFrame, 0, _frameSize>>2);
  memset(quantVFrame, 0, _frameSize>>2);

  _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
  _trans->quantization(dctUFrame, quantUFrame, cw, ch);
  _trans->quantization(dctVFrame, quantVFrame, cw, ch);
  _numChnCodeBands = 0;
  bitsU_in = _cavlc->encode(quantUFrame, _bsU);
  bitsV_in = _cavlc->encode(quantVFrame, _bsV);
# if HARDWARE_FLOW
  if (bitsU_in%32 != 0) {
    int dummy = 32 - (bitsU_in%32);
    _bsU->write(0, dummy);
  }
  if (bitsV_in%32 != 0) {
    int dummy = 32 - (bitsV_in%32);
    _bsV->write(0, dummy);
  }
# endif // HARDWARE_FLOW
  cout << "Bits written (U): " << bitsU_in << endl;
  cout << "Bits written (V): " << bitsV_in << endl;
  _bsU->flush();
  _bsV->flush();
  _files->getFile("wzU_in")->closeFile();
  _files->getFile("wzV_in")->closeFile();
}

void EncTest::cavlcTest1() {
  cout << "####################################" << endl; 
  cout << "FILE->cavlcDec:" << endl;
  cout << "####################################" << endl; 
  cout << "Quantized vs. cavlcDecoded:" << endl;
  memset(dctUFrame, 0, _frameSize>>2);
  memset(dctVFrame, 0, _frameSize>>2);
  memset(quantUFrame, 0, _frameSize>>2);
  memset(quantVFrame, 0, _frameSize>>2);
  memset(iDecodedU, 0, _frameSize>>2);
  memset(iDecodedV, 0, _frameSize>>2);
  string ubs = "test.u.bin";
  string vbs = "test.v.bin";
  _files->addFile("wzU_out", ubs.c_str())->openFile("rb");
  _files->addFile("wzV_out", vbs.c_str())->openFile("rb");
  Bitstream* bsU_out = new Bitstream(1024, _files->getFile("wzU_out")->getFileHandle());
  Bitstream* bsV_out = new Bitstream(1024, _files->getFile("wzV_out")->getFileHandle());
  CavlcDec* cavlcDecU = new CavlcDec(this, 4);
  CavlcDec* cavlcDecV = new CavlcDec(this, 4);

  _trans->dctTransform(chromaResidue, dctUFrame, cw, ch);
  _trans->dctTransform(chromaResidue + chsize, dctVFrame, cw, ch);
  _trans->quantization(dctUFrame, quantUFrame, cw, ch);
  _trans->quantization(dctVFrame, quantVFrame, cw, ch);

  int bitsU_out = 0;
  int bitsV_out = 0;
  int buff1[16];
  int buff2[16];
  // read bits from bitstream
  for (int j = 0; j < ch; j += 4) {
    for (int i = 0; i < cw; i += 4) {
      bitsU_out += cavlcDecU->decode(iDecodedU, i, j, bsU_out);
      bitsV_out += cavlcDecV->decode(iDecodedV, i, j, bsV_out);

      for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
          buff1[y*4+x] = quantVFrame[(j+y)*cw+i+x];
          buff2[y*4+x] = iDecodedV[(j+y)*cw+i+x];
        }
      }
      assert(calcMSE(buff1, buff2, 16) == 0);
    }
  }
  cout << "bits read (U):" << bitsU_out << endl;
  cout << "bits read (V):" << bitsV_out << endl;
  assert(bitsU_in == bitsU_out);
  assert(bitsV_in == bitsV_out);
  assert(calcMSE(quantUFrame, iDecodedU, chsize) == 0);
  cout << (calcMSE(quantVFrame, iDecodedV, chsize)) << endl;//== 0);
  cout << (calcPSNR(quantVFrame, iDecodedV, chsize)) << endl;//== 0);
}

//void EncTest::cavlcTest2() {
//  cout << "####################################" << endl; 
//  cout << "FILE->cavlcDec:" << endl;
//  cout << "####################################" << endl; 
//  cout << "Quantized vs. cavlcDecoded:" << endl;
//  memset(iDecodedU, 0, _frameSize>>2);
//  memset(iDecodedV, 0, _frameSize>>2);
//  memset(iDctU, 0, _frameSize>>2);
//  memset(iDctV, 0, _frameSize>>2);
//  memset(iQuantU, 0, _frameSize>>2);
//  memset(iQuantV, 0, _frameSize>>2);
//  string ubs = "test.u.bin";
//  string vbs = "test.v.bin";
//  _files->addFile("wzU_out", ubs.c_str())->openFile("rb");
//  _files->addFile("wzV_out", vbs.c_str())->openFile("rb");
//  Bitstream* bsU_out = new Bitstream(1024, _files->getFile("wzU_out")->getFileHandle());
//  Bitstream* bsV_out = new Bitstream(1024, _files->getFile("wzV_out")->getFileHandle());
//  CavlcDec* cavlcDecU = new CavlcDec(this, 4);
//  CavlcDec* cavlcDecV = new CavlcDec(this, 4);
//
//  // cavlc decode
//  int bitsU_out = 0;
//  int bitsV_out = 0;
//  for (int j = 0; j < ch; j += 4) {
//    for (int i = 0; i < cw; i += 4) {
//      bitsU_out += cavlcDecU->decode(iDecodedU, i, j, bsU_out);
//      bitsV_out += cavlcDecV->decode(iDecodedV, i, j, bsV_out);
//    }
//  }
//
//  // inverse Q & DCT
//  _trans->invQuantization(iDecodedU, iQuantU, cw, ch);
//  _trans->invQuantization(iDecodedV, iQuantV, cw, ch);
//  _trans->invDctTransform(iQuantU, iDctU, cw, ch);
//  _trans->invDctTransform(iQuantV, iDctV, cw, ch);
//
//  double currPSNR = calcPSNR(currChroma, iQuantU, chsize);
//  cout << "PSNR (U): " << currPSNR << endl;
//  currPSNR = calcPSNR(currChroma, iQuantV, chsize);
//  cout << "PSNR (V): " << currPSNR << endl;
//
//  // visualize
//  int* container[_frameSize>>1];
//  memcpy(container, quantUFrame, chsize);
//  memcpy(container+chsize, iDecodedU, chsize);
//  cv::Mat mu(ch*2, cw, CV_8UC1, container);
//
////  double alpha = 1.0;
////  int beta = 20;
////  for (int j = 0; j < 2*ch; j++) {
////    for (int i = 0; i < cw; i++) {
////      mu.at<int>(j,i) = cv::saturate_cast<uchar>(alpha*mu.at<int>(j,i) + beta);
////    }
////  }
//
//  cv::imshow("Frame 1", mu);
//  cv::waitKey(0);
//}



// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc != 2) {
    cerr << endl;
    cerr << "Usage: ./tester ";
    cerr << "[config file]";
    cerr << endl;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1]);
    EncTest* test = new EncTest(configMap);

    test->transTest();
    test->qTest1();
    test->qTest2();
    test->cavlcEnc();
    test->cavlcTest1();

  }
  return 0;
}
