#ifndef calculations_h_
#define calculations_h_

#include <cmath>

#include "defs.h"

int calcSAD(imgpel* blk1, imgpel* blk2, int width, int blocksize);

int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2,
            int s1,int s2,int blocksize);

int calcSAD(imgpel* blk1, imgpel* blk2, int px, int py, 
            int rx,int ry, const int blocksize, int width, int height);

int calcDist(imgpel* blk1, imgpel* blk2, int width1, int width2,
             int s1, int s2, int blocksize);


// Templated functions
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T, class U>
double calcMSE(T* img1, U* img2,int length)
{
  float MSE=0;
  for(int i=0;i<length;i++)
    {
      MSE+=pow(float(img1[i]-img2[i]),float(2.0))/length;
    }
  return MSE;
}
template<class T, class U>
double calcPSNR(T* img1, U* img2,int length)
{
  return 10*log10(255*255/calcMSE(img1, img2, length));
}
#endif
