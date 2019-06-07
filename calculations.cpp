#include <cstdlib>

#include "calculations.h"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcDist(imgpel* blk1, imgpel* blk2, int width1, int width2,
             int s1, int s2, int blocksize){
  int iDist=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+s1*x+s1*y*width1);
      imgpel pel2=*(blk2+s2*x+s2*y*width2);
      iDist+=(pel1-pel2)*(pel1-pel2);
    }
  return iDist;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcSAD(imgpel* blk1, imgpel* blk2, int width, int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+x+y*width);
      imgpel pel2=*(blk2+x+y*width);
      sad+=abs(pel1-pel2);
    }
  return sad;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2, 
            int s1,int s2,int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+s1*x+s1*y*width1);
      imgpel pel2=*(blk2+s2*x+s2*y*width2);
      sad+=abs(pel1-pel2);
    }
  return sad;
}
