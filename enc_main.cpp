
#include <iostream>
#include <exception>

#include "codec.h"
#include "encoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 3) {
    cerr << endl;
    cerr << "Usage: ./enDVC ";
    cerr << "[config file] ";
    cerr << "[input video]";
    cerr << endl;
    return 1;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1], true);
    configMap["SrcFile"] = argv[2];
    if (FILE *file = fopen(configMap["SrcFile"].c_str(), "r"))
      fclose(file);
    else
      throw invalid_argument("Invalid source file");
    Encoder* encoder = new Encoder(configMap);

    encoder->encodeKeyFrame();
    encoder->encodeWzHeader();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  return 0;
}

