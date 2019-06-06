
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "transform.h"
#include "codec.h"
#include "fileManager.h"

const int Transform::DctScaler[4][4] = {
  { 1, 81,  1, 81},
  {81, 51, 81, 51},
  { 1, 81,  1, 81},
  {81, 51, 81, 51}
};

const int Transform::DctShift[4][4] = {
  {2, 9, 2, 9},
  {9, 9, 9, 9},
  {2, 9, 2, 9},
  {9, 9, 9, 9}
};

const int Transform::IdctScaler[4][4] = {
  {64, 81, 64, 81},
  {81, 51, 81, 51},
  {64, 81, 64, 81},
  {81, 51, 81, 51}
};

const int Transform::IdctShift[4][4] = {
  {8, 8, 8, 8},
  {8, 7, 8, 7},
  {8, 8, 8, 8},
  {8, 7, 8, 7}
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template void Transform::dctTransform(int*    src, int* dst, bool isChr);
template void Transform::dctTransform(imgpel* src, int* dst, bool isChr);

template <typename T, typename U>
void Transform::dctTransform(T* src, U* dst, bool isChr)
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }
  int** srcBuffer = new int*[4];
  int** dstBuffer = new int*[4];
  for (int i = 0; i < 4; i++) {
    srcBuffer[i] = new int[4];
    dstBuffer[i] = new int[4];
  }
  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4) {
# if INTEGER_DCT
      //copy
      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++)
          srcBuffer[j][i] = (int)src[(x+i) + (y+j)*width];

      forward4x4(srcBuffer, dstBuffer, 0, 0);

      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++) {
          int sign = (dstBuffer[j][i] >= 0) ? 1 : -1;

          dst[(x+i)+(y+j)*width] = sign*((abs(dstBuffer[j][i])*DctScaler[j][i] +
                                         (0x1<<(DctShift[j][i]-1))) >> DctShift[j][i]);
        }
# else // if !INTEGER_DCT
      dct4x4(src, dst, x, y);
# endif // INTEGER_DCT
    }
  for (int i = 0; i < 4; i++) {
    delete [] srcBuffer[i];
    delete [] dstBuffer[i];
  }
  delete [] srcBuffer;
  delete [] dstBuffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::dctTransform(float* src, float* dst, bool isChr)
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }
  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4)
      dct4x4(src, dst, x, y, width);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template void Transform::dct4x4(float*  src, float* dst, int, int, int);
template void Transform::dct4x4(int*    src, int*   dst, int, int, int);
template void Transform::dct4x4(imgpel* src, int*   dst, int, int, int);

template <typename T, typename U>
void Transform::dct4x4(T* src, U* dst, int x, int y, int width)
{
  float   c;
  float buffer[4][4];

  for (int v = y; v < y+4; v++) {
    if (v == y)
      c = 1/2.0;
    else
      c = (float)(sqrtf(2.0)/2.0);

    for (int i = x; i < x+4; i++) {
      float t = 0.0;

      for (int j = y; j < y+4 ;j++)
        t += c*(float)cos((2*(j-y)+1)*(v-y)*PI/8.0)*(float)src[i+j*width];

      buffer[v-y][i-x] = t;
    }
  }

  for (int u = x; u < x+4; u++) {
    if (u == x)
      c = 1/2.0;
    else
      c = (float)(sqrtf(2.0)/2.0);

    for (int v = y; v < y+4; v++) {
      float t = 0.0;

      for (int i = x; i < x+4; i++)
        t += c*(float)cos((2*(i-x)+1)*(u-x)*PI/8.0)*buffer[v-y][i-x];

      dst[u+v*width] = (U)t;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::forward4x4(int** src, int** dst, int x, int y)
{
  int   buffer[16];
  int*  ptr = buffer;

  int   p0, p1, p2, p3;
  int   t0, t1, t2, t3;

  // Horizontal
  for (int i = y; i < y + BLOCK_SIZE; i++) {
    p0 = src[i][x + 0];
    p1 = src[i][x + 1];
    p2 = src[i][x + 2];
    p3 = src[i][x + 3];

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    *(ptr++) =  t0 + t1;
    *(ptr++) = (t3 << 1) + t2;
    *(ptr++) =  t0 - t1;
    *(ptr++) =  t3 - (t2 << 1);
  }

  // Vertical
  for (int i = 0; i < BLOCK_SIZE; i++) {
    ptr = buffer + i;

    p0 = *ptr;
    p1 = *(ptr += BLOCK_SIZE);
    p2 = *(ptr += BLOCK_SIZE);
    p3 = *(ptr += BLOCK_SIZE);

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    dst[y + 0][x + i] = t0 +  t1;
    dst[y + 1][x + i] = t2 + (t3 << 1);
    dst[y + 2][x + i] = t0 -  t1;
    dst[y + 3][x + i] = t3 - (t2 << 1);
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template void Transform::invDctTransform(int* src, int* dst, bool isChr);
template void Transform::invDctTransform(int* src, imgpel* dst, bool isChr);

template <typename T>
void Transform::invDctTransform(int* src, T* dst, bool isChr)
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }
  int** srcBuffer = new int*[4];
  int** dstBuffer = new int*[4];

  for (int i = 0; i < 4; i++) {
    srcBuffer[i] = new int[4];
    dstBuffer[i] = new int[4];
  }
  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4) {
# if INTEGER_DCT

      //copy
      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++) {
          int sign = (src[(x+i) + (y+j)*width] >= 0)? 1 : -1;

          srcBuffer[j][i] = sign*((abs(src[x+i+(y+j)*width])*IdctScaler[j][i]+
                                   (0x1<<(IdctShift[j][i]-7))) >> (IdctShift[j][i]-6));
        }

      inverse4x4(srcBuffer, dstBuffer, 0, 0);

      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++) {
          int level;

          if (dstBuffer[j][i] >= 0)
            level = (dstBuffer[j][i]+32)>>6;
          else
            level = (dstBuffer[j][i]-32)>>6;

          dst[(x+i) + (y+j)*width] = level;
        }
# else // if !INTEGER_DCT
      idct4x4(src, dst, x, y);
# endif // INTEGER_DCT
    }
  for (int i = 0; i < 4; i++) {
    delete [] srcBuffer[i];
    delete [] dstBuffer[i];
  }
  delete [] srcBuffer;
  delete [] dstBuffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::inverse4x4(int** src, int** dst, int x, int y)
{
  int   buffer[16];
  int*  ptr = buffer;

  int   p0, p1, p2, p3;
  int   t0, t1, t2, t3;

  // Horizontal
  for (int i = y; i < y + BLOCK_SIZE; i++) {
    t0 = src[i][x + 0];
    t1 = src[i][x + 1];
    t2 = src[i][x + 2];
    t3 = src[i][x + 3];

    p0 =  t0 + t2;
    p1 =  t0 - t2;
    p2 = (t1 >> 1) - t3;
    p3 =  t1 + (t3 >> 1);

    *(ptr++) = p0 + p3;
    *(ptr++) = p1 + p2;
    *(ptr++) = p1 - p2;
    *(ptr++) = p0 - p3;
  }

  //  Vertical
  for (int i = 0; i < BLOCK_SIZE; i++) {
    ptr = buffer + i;

    t0 = *ptr;
    t1 = *(ptr += BLOCK_SIZE);
    t2 = *(ptr += BLOCK_SIZE);
    t3 = *(ptr += BLOCK_SIZE);

    p0 =  t0 + t2;
    p1 =  t0 - t2;
    p2 = (t1 >> 1) - t3;
    p3 =  t1 + (t3 >> 1);

    dst[y    ][x + i] = p0 + p3;
    dst[y + 1][x + i] = p1 + p2;
    dst[y + 2][x + i] = p1 - p2;
    dst[y + 3][x + i] = p0 - p3;
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template void Transform::idct4x4(int* src, int*    dst, int x, int y, int);
template void Transform::idct4x4(int* src, imgpel* dst, int x, int y, int);

template <typename T>
void Transform::idct4x4(int* src, T* dst, int x, int y, int width)
{
  float   c;
  float buffer[4][4];

  for (int v = y; v < y+4; v++) {
    for (int i = x; i < x+4; i++) {
      float t = 0.0;

      for (int j = y; j < y+4; j++) {
        if (j == y)
          c = 1/2.0;
        else
          c = (float)(sqrtf(2.0)/2.0);

        t += c*(float)cos((2*(v-y)+1)*(j-y)*PI/8.0)*(float)src[i+j*width];
      }

      buffer[v-y][i-x] = t;
    }
  }

  for (int u = x; u < x+4; u++) {
    for (int v = y; v < y+4; v++) {
      float t = 0.0;

      for (int i = x; i < x+4; i++) {
        if (i == x)
          c = 1/2.0;
        else
          c = (float)(sqrtf(2.0)/2.0);

        t += c*(float)cos((2*(u-x)+1)*(i-x)*PI/8.0)*buffer[v-y][i-x];
      }

      dst[u+v*width] = (T)t;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::quantization(int* src, int* dst, bool isChr)
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }
  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4)
      quan4x4(src, dst, x, y, isChr);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::quan4x4(int* src, int* dst, int x, int y, bool isChr)
{
  int index;
  int qp, width;
  if (isChr) {
    qp = _codec->getChrQp();
    width = _codec->getFrameWidth()>>1;
  } else {
    qp = _codec->getQp();
    width = _codec->getFrameWidth();
  }

  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++) {
      index = (x+i) + (y+j)*width;

      if (i == 0 && j == 0) {
        dst[index] = (abs(src[index])/_codec->getQuantStep(i, j, isChr));

        if (src[index]/(_codec->getQuantStep(i, j, isChr)) < 0)
          dst[index] |= (0x1 << (_codec->getQuantMatrix(qp, i, j)-1));
      }
      else if (_codec->getQuantMatrix(qp, i, j) != 0) {
        dst[index] = (abs(src[index])/_codec->getQuantStep(i, j, isChr));

        if ((src[index]/_codec->getQuantStep(i, j, isChr)) < 0)
          dst[index] |= (0x1 << (_codec->getQuantMatrix(qp, i, j)-1));
      }
      else {
        dst[index] = 0;
        break;
      }
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
#if SI_REFINEMENT
void Transform::invQuantization(int* src, int* dst, int* si,
                                int offsetX, int offsetY, bool isChr)
#else
void Transform::invQuantization(int* src, int* dst, int* si, bool isChr)
#endif
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }

  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4) {
#if SI_REFINEMENT
      invquan4x4(src, dst, si, x, y, offsetX, offsetY, isChr);
#else
      invquan4x4(src, dst, si, x, y, isChr);
#endif
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
#if SI_REFINEMENT
void Transform::invQuantization(int* src, int* dst,
                                int offsetX, int offsetY, bool isChr)
#else
void Transform::invQuantization(int* src, int* dst, bool isChr)
#endif
{
  int width, height;
  if (isChr) {
    width = _codec->getFrameWidth()>>1;
    height = _codec->getFrameHeight()>>1;
  } else {
    width = _codec->getFrameWidth();
    height = _codec->getFrameHeight();
  }
  for (int y = 0; y < height; y += 4)
    for (int x = 0; x < width; x += 4) {
#if SI_REFINEMENT
      invquan4x4(src, dst, x, y, offsetX, offsetY, isChr);
#else
      invquan4x4(src, dst, x, y, isChr);
#endif
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
# if SI_REFINEMENT
void Transform::invquan4x4(int* src, int* dst, int* si, int x, int y,
                           int i, int j, bool isChr)
# else
void Transform::invquan4x4(int* src, int* dst, int* si, int x, int y, bool isChr)
# endif
{
  int qp, width;
  if (isChr) {
    qp = _codec->getChrQp();
    width = _codec->getFrameWidth()>>1;
  } else {
    qp = _codec->getQp();
    width = _codec->getFrameWidth();
  }

# if !SI_REFINEMENT
  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++)
# endif
    {
      int     index = (x+i) + (y+j)*width;
      int     step  = _codec->getQuantStep(i, j, isChr);
      int     mask  = (0x1<<(_codec->getQuantMatrix(qp, i, j)-1)) - 1;
      int     sign  = (src[index]>>(_codec->getQuantMatrix(qp, i, j)-1)) & 0x1;
      int     value = src[index] & mask;
      double  alpha = _codec->getAlpha()[index];

      if (i == 0 && j == 0) {
        dst[index] = (sign == 0) ? value*step : -value*step;

        mmse(dst+index, si+index, dst+index, step, alpha);
      }
      else if (_codec->getQuantMatrix(qp, i, j) != 0) {
        dst[index] = (sign == 0) ? value*step : -value*step;

        mmse(dst+index, si+index, dst+index, step, alpha);
      }
      else
        dst[index] = si[index];
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
# if SI_REFINEMENT
void Transform::invquan4x4(int* src, int* dst, int x, int y,
                           int i, int j, bool isChr)
# else
void Transform::invquan4x4(int* src, int* dst, int x, int y, bool isChr)
# endif
{
  int qp, width;
  if (isChr) {
    qp = _codec->getChrQp();
    width = _codec->getFrameWidth()>>1;
  } else {
    qp = _codec->getQp();
    width = _codec->getFrameWidth();
  }

# if !SI_REFINEMENT
  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++)
# endif
    {
      int index = (x+i) + (y+j)*width;
      int step  = _codec->getQuantStep(i, j, isChr);
      int mask  = (0x1<<(_codec->getQuantMatrix(qp, i, j)-1)) - 1;
      int sign  = (src[index]>>(_codec->getQuantMatrix(qp, i, j)-1)) & 0x1;
      int value = src[index] & mask;

      dst[index] = (sign == 0) ? value*step : -value*step;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Transform::mmse(int* x, int* y, int* z, int step, double alpha)
{
  double gamma;
  double delta;

  if (*z > 0) {
    gamma = *y - *z;
    delta = *z+step - *y;

    if (*y < *z)
      *x = (int)((double)*z + 1.0/alpha + (double)step/(1-exp(alpha*step)));
    else if (*y >= *z+step)
      *x = (int)((double)(*z + step) - 1.0/alpha - (double)step/(1-exp(alpha*step)));
    else
      *x = (int)((double)*y + ((gamma+1.0/alpha)*exp(-alpha*gamma) - (delta+1.0/alpha)*exp(-alpha*delta)) /
                              (2 - (exp(-alpha*gamma) + exp(-alpha*delta))));
  }
  else if (*z < 0) {
    gamma = *z - *y;
    delta = *y - (*z-step);

    if (*y >= *z)
      *x = (int)((double)*z - 1.0/alpha - (double)step/(1-exp(alpha*step)));
    else if (*y <= *z-step)
      *x = (int)((double)(*z - step) + 1.0/alpha + (double)step/(1-exp(alpha*step)));
    else
      *x = (int)((double)*y - ((gamma+1.0/alpha)*exp(-alpha*gamma) - (delta+1.0/alpha)*exp(-alpha*delta)) /
                              (2 - (exp(-alpha*gamma) + exp(-alpha*delta))));
  }
  else {
    if (*y >= step)
      *x = (int)((double)step - 1.0/alpha - (double)2*step/(1-exp(alpha*2*step)));
    else if (*y <= -step)
      *x = (int)(-step + 1.0/alpha + (double)2*step/(1-exp(alpha*2*step)));
    else
      *x = *y;
  }
}

