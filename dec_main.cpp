
#include <iostream>

#include "decoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 2) {
    cerr << "Usage: ./deDVC [config file]" << endl;
    return 1;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1]);
    Decoder* decoder = new Decoder(configMap);

    decoder->decodeWZframe();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  //system("pause");
  return 0;
}

