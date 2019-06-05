
#include <iostream>

#include "decoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 3) {
    cerr << "Usage: ./deDVC [config file] [input video]" << endl;
    return 1;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1], false);
    configMap["SrcFile"] = argv[2];
    if (FILE *file = fopen(configMap["SrcFile"].c_str(), "r"))
      fclose(file);
    else
      throw invalid_argument("Invalid source file");
    Decoder* decoder = new Decoder(configMap);

    decoder->decodeWzFrame();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  return 0;
}

