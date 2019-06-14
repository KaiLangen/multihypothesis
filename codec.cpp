
#include <sstream>
#include <iostream>
#include <fstream>
#include <exception>

#include "codec.h"

const int Codec::ResidualBlockSize = 8;
const int Codec::SkipBlockSize = 4;

const int Codec::QuantMatrix[8][4][4] = {
 {{4, 3, 0, 0},
  {3, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0}},
 {{5, 3, 0, 0},
  {3, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0}},
 {{5, 3, 2, 0},
  {3, 2, 0, 0},
  {2, 0, 0, 0},
  {0, 0, 0, 0}},
 {{5, 4, 3, 2},
  {4, 3, 2, 0},
  {3, 2, 0, 0},
  {2, 0, 0, 0}},
 {{5, 4, 3, 2},
  {4, 3, 2, 2},
  {3, 2, 2, 0},
  {2, 2, 0, 0}},
 {{6, 4, 3, 3},
  {4, 3, 3, 2},
  {3, 3, 2, 0},
  {3, 2, 0, 0}},
 {{6, 5, 4, 3},
  {5, 4, 3, 2},
  {4, 3, 2, 0},
  {3, 2, 0, 0}},
 {{7, 6, 5, 4},
  {6, 5, 4, 3},
  {5, 4, 3, 0},
  {4, 3, 0, 0}}
};

const int Codec::BitPlaneNum[8] = {10, 11, 17, 30, 36, 41, 50, 59};

const int Codec::MaxBitPlane[4][4] = {
  {10, 9, 8, 8},
  { 9, 8, 8, 7},
  { 8, 8, 7, 6},
  { 8, 7, 6, 6}
};

const int Codec::MinQStepSize[8][4][4] = {
 {{32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32}},
 {{32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32}},
 {{32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32}},
 {{32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32},
  {32, 32, 32, 32}},
 {{16, 16, 16, 16},
  {16, 16, 16, 16},
  {16, 16, 16, 16},
  {16, 16, 16, 16}},
 {{16, 16, 16, 16},
  {16, 16, 16, 16},
  {16, 16, 16, 16},
  {16, 16, 16, 16}},
 {{10, 10, 10, 12},
  {10, 10, 12, 16},
  {10, 12, 16, 16},
  {12, 16, 16, 16}},
 {{ 8,  8,  8, 10},
  { 8,  8, 10, 16},
  { 8, 10, 16, 18},
  {10, 16, 18, 18}}
};

// Zigzag scan order
const int Codec::ScanOrder[16][2] = {
  {0, 0},
  {1, 0},
  {0, 1},
  {0, 2},
  {1, 1},
  {2, 0},
  {3, 0},
  {2, 1},
  {1, 2},
  {0, 3},
  {1, 3},
  {2, 2},
  {3, 1},
  {3, 2},
  {2, 3},
  {3, 3}
};

const int Codec::HuffmanCodeValue[4][3][16] = {
  // QP = 1 or 2
  {{  3,  1,  4,  0,  3, 23, 21,  4, 40, 11, 83, 82, 21, 41, 40, 22},   // type 0
   {  3,  0,  3,  9,  4, 23, 16, 11, 45, 35, 21, 34, 20, 89, 88, 10},   // type 1
   {  2,  2,  3, 15, 13, 12,  5,  4,  3,  2,  0, 28, 29,  3,  2,  3}},  // type 2
  // QP = 3 or 4
  {{  3,  1,  4, 11,  3,  1, 20,  4,  0, 11,  3,  2, 20, 43, 42, 21},   // type 0
   {  3,  0,  3,  9, 23, 21, 17, 45, 41, 33, 89, 32, 88, 81, 80,  2},   // type 1
   {  2,  3,  5,  2,  0,  3,  6,  2, 19, 15, 14, 17, 16, 37, 36,  3}},  // type 2
  // QP = 5 or 6
  {{  3,  1,  4, 11,  3,  1, 20,  1, 11, 10,  8,  1,  0, 19, 18, 21},   // type 0
   {  0,  7, 13, 10, 25, 23, 17, 49, 45, 33, 32, 97, 96, 89, 88,  9},   // type 1
   {  2, 14, 31, 25, 60, 53, 54, 49,122, 48,123,111,110,105,104,  0}},  // type 2
  // QP = 7 or 8
  {{  0,  7,  4, 12, 27, 23, 20, 52, 44,106, 91,215,214,181,180, 21},   // type 0
   {  0,  6, 15,  9, 29, 21, 16, 57, 41, 35, 34,113,112, 81, 80, 11},   // type 1
   {  0,  4, 11, 30, 21, 20, 62, 56, 58,127,126,119,115,114,118,  6}}   // type 2
};

const int Codec::HuffmanCodeLength[4][3][16] = {
  // QP = 1 or 2
  {{  2,  2,  3,  3,  4,  5,  5,  5,  6,  6,  7,  7,  7,  8,  8,  5},   // type 0
   {  2,  2,  3,  4,  4,  5,  5,  5,  6,  6,  6,  6,  6,  7,  7,  4},   // type 1
   {  2,  3,  4,  5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  2}},  // type 2
  // QP = 3 or 4
  {{  2,  2,  3,  4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  8,  8,  5},   // type 0
   {  2,  2,  3,  4,  5,  5,  5,  6,  6,  6,  7,  6,  7,  7,  7,  3},   // type 1
   {  2,  3,  4,  4,  4,  5,  5,  5,  6,  6,  6,  6,  6,  7,  7,  2}},  // type 2
  // QP = 5 or 6
  {{  2,  2,  3,  4,  4,  4,  5,  5,  6,  6,  6,  6,  6,  7,  7,  5},   // type 0
   {  1,  3,  4,  4,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  4},   // type 1
   {  2,  4,  5,  5,  6,  6,  6,  6,  7,  6,  7,  7,  7,  7,  7,  1}},  // type 2
  // QP = 7 or 8
  {{  1,  3,  3,  4,  5,  5,  5,  6,  6,  7,  7,  8,  8,  8,  8,  5},   // type 0
   {  1,  3,  4,  4,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  4},   // type 1
   {  1,  3,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  7,  7,  7,  3}}   // type 2
};

void Codec::getRecFrame(imgpel* recon, imgpel* bRef, imgpel* fRef,
                        int* curr,  int* rcList, bool isChr)
{
  int blockCount = 0;
  int width, height;
  if (isChr) {
    width = _frameWidth>>1;
    height = _frameHeight>>1;
  } else {
    width = _frameWidth;
    height = _frameHeight;
  }
  imgpel* refFrame;

  for (int j = 0; j < height; j += ResidualBlockSize)
    for (int i = 0; i < width; i += ResidualBlockSize) {
# if HARDWARE_OPT
      refFrame = bRef;
      (void) fRef;
      (void) rcList;
# else // if !HARDWARE_OPT
      refFrame = (rcList[blockCount] == 0) ? bRef : fRef;
#endif
      for (int y = 0; y < ResidualBlockSize; y++)
        for (int x = 0; x < ResidualBlockSize; x++) {
          int idx = (i+x) + (j+y)*width;

          recon[idx] = curr[idx] + refFrame[idx];
        }
      blockCount++;
    }
}

using namespace std;

map<string, string>&
readConfig(string filename, bool isEnc)
{
  string line;
  ifstream cfile(filename);
  static map<string, string> configMap;
  if (cfile)
  {
    while( getline(cfile, line) )
    {
      string key, value;
      istringstream is_line(line);

      // Check if the first non-whitespace is a #
      if( getline(is_line >> ws, key, '=')  && !key.empty() && key[0] != '#')
        if( getline(is_line >> ws, value) )
          configMap[key] = value;
    }

    // check the validity of inputs
    if (isEnc) { // encoder specific inputs
      int qp = atoi(configMap["WzQP"].c_str());
      if (qp < 0 || qp >= 8)
        throw invalid_argument("Invalid quantization paramater");
      int chrQp = atoi(configMap["ChrQP"].c_str());
      if (chrQp < 0 || chrQp >= 8)
        throw invalid_argument("Invalid Chroma QP");
      int keyQp = atoi(configMap["KeyQP"].c_str());
      if (keyQp < 0 || keyQp >= 52)
        throw invalid_argument("Invalid Key-frame QP");
      int numFrames = atoi(configMap["NumFrames"].c_str());
      if (numFrames <= 0)
        throw invalid_argument("Invalid number of frames");
      int gop = atoi(configMap["Gop"].c_str());
      if (gop < 2 || gop % 2 == 1)
        throw invalid_argument("Invalid GOP size");
      string seqType = configMap["SequenceType"];
      if ((seqType.compare("CIF") != 0) && (seqType.compare("QCIF") != 0))
        throw invalid_argument("Invalid GOP size");
      if (FILE *file = fopen(configMap["KeyFile"].c_str(), "w"))
        fclose(file);
      else
        throw invalid_argument("Unable to create key-frame file");

    } else { // decoder specific inputs
      int blockSize = atoi(configMap["BlockSize"].c_str());
      if (blockSize < 8 || blockSize > 32)
        throw invalid_argument("Invalid search-block size");
      int searchWindow = atoi(configMap["SearchWindowSize"].c_str());
      if (searchWindow < 0 || searchWindow > 32)
        throw invalid_argument("Invalid search-window size");
      int nRefFrames = atoi(configMap["NumRefFrames"].c_str());
      if (nRefFrames < 0)
        throw invalid_argument("Invalid frame sliding-window size");
      int spatialSmoothing = atoi(configMap["SearchWindowSize"].c_str());
      if (spatialSmoothing < 0)
        throw invalid_argument("Invalid spatial smoothing parameter");
      int meMode = atoi(configMap["MEMode"].c_str());
      if (meMode < 0 || meMode > 1)
        throw invalid_argument("Invalid Motion Estimation Mode: must use 0 for Chroma-ME or 1 for oracle");
      if (FILE *file = fopen(configMap["WZFile"].c_str(), "r"))
        fclose(file);
      else
        throw invalid_argument("Unable to read WZ binary file");
      if (FILE *file = fopen(configMap["KeyFile"].c_str(), "r"))
        fclose(file);
      else
        throw invalid_argument("Unable to read key-frame binary file");
    }

    cfile.close();
    return configMap;
  }
  else
  {
    cerr << "No such file: " << filename << endl;
    throw invalid_argument("Invalid config file");
  }
}
