#include <map>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "config.h"
#include "calculations.h"
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
  _chrQp       = atoi(configMap["ChrQP"].c_str()); 
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
  _files->addFile("oracle", configMap["OracleFile"]);

  string ubs = wzFileName.substr(0, wzFileName.find(".bin")) + ".u.bin";
  string vbs = wzFileName.substr(0, wzFileName.find(".bin")) + ".v.bin";
  _files->addFile("wzU",  ubs.c_str())->openFile("wb");
  _files->addFile("wzV",  vbs.c_str())->openFile("wb");

  _bs = new Bitstream(1024, _files->getFile("wz")->getFileHandle());
  _bsU = new Bitstream(1024, _files->getFile("wzU")->getFileHandle());
  _bsV = new Bitstream(1024, _files->getFile("wzV")->getFileHandle());
  _frameSize        = _frameWidth * _frameHeight;
  _bitPlaneLength   = _frameSize / 16;
  _fb = new FrameBuffer(_frameWidth, _frameHeight);
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
  cmd << "-p FramesToBeEncoded=" << ((_numFrames/_gop)+1) << " ";
  cmd << "-p QPISlice=" << _keyQp << " ";
  cmd << "-p FrameSkip=" << _gop-1 << " ";
  cmd << "-p SourceWidth=" << _frameWidth << " ";
  cmd << "-p SourceHeight=" << _frameHeight << " ";
  cmd << " > jm.log;";
  cmd << "cp stats.dat " << BIN_DIR;
  system(cmd.str().c_str());
  _files->getFile("key")->openFile("rb");

  cout << "Running JM to encode oracle frames" << endl;
  string oracleFileName = _files->getFile("oracle")->getFileName();
  stringstream cmd2(stringstream::in | stringstream::out);
  cmd2 << "cd jm; ";
  cmd2 << "./lencod.exe -d encoder_intra_main.cfg ";
  cmd2 << "-p InputFile=\"" << BIN_DIR << "/" << srcFileName << "\" ";
  cmd2 << "-p ReconFile=\"" << BIN_DIR << "/" << oracleFileName << "\" ";
  cmd2 << "-p FramesToBeEncoded=" << _numFrames << " ";
  cmd2 << "-p QPISlice=" << _keyQp << " ";
  cmd2 << "-p SourceWidth=" << _frameWidth << " ";
  cmd2 << "-p SourceHeight=" << _frameHeight << " ";
  cmd2 << " > jm.log;";
  cmd2 << "cp stats.dat " << BIN_DIR;
  system(cmd2.str().c_str());
  _files->getFile("oracle")->openFile("rb");
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeWzHeader()
{
  // Encode width and height in terms of macroblock
  _bs->write(_frameWidth/16, 8);
  _bs->write(_frameHeight/16, 8);
  _bs->write(_qp, 8);
  _bs->write(_chrQp, 8);
  _bs->write(_numFrames, 16);
  _bs->write(_gop, 8);
  _bs->flush();
}

