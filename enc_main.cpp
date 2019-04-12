
#include <iostream>

#include "encoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 2) {
    cerr << endl;
    cerr << "Usage: ./enDVC ";
    cerr << "[config file]";
    cerr << endl;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1]);
    Encoder* encoder = new Encoder(configMap);

    encoder->encodeKeyFrame();
    encoder->encodeWzFrame();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  //system("pause");
  return 0;
}

