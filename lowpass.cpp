#include "types.h"
#include "obme.h"

void lowpassFilter(imgpel* src, imgpel* dst,
                   const int boxSize, const int width, const int height)
{
  const int left = boxSize/2;
  const int right = boxSize - (left + 1);

  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++) {
      int sum   = 0;
      int total = 0;
  
      // Sum up pixels within the box
      for (int j = y-left; j <= y+right; j++)
        for (int i = x-left; i <= x+right; i++)
          if (i >= 0 && i < width &&
              j >= 0 && j < height) {
            sum += src[i+j*width];
            total++;
          }
 
      dst[x+y*width] = (imgpel)(sum/total);
    }
}
