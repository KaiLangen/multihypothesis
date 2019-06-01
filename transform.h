
#ifndef ENCODER_INC_TRANSFORM_H
#define ENCODER_INC_TRANSFORM_H

#include "defs.h"

class Codec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Transform
{
public:
  Transform(Codec* codec) { _codec = codec;}

  // DCT
  template <typename T, typename U>
  void dctTransform(T* src, U* dst, int w, int h);

  void dctTransform(float* src, float* dst, int w, int h);

  template <typename T, typename U>
  void dct4x4(T* src, U* dst, int x, int y, int w);

  // Inverse DCT
  template <typename T>
  void invDctTransform(int* src, T* dst, int w, int h);

  // Quantization
  void quantization(int* src, int* dst, int w, int h);

  // Inverse quantization
# if SI_REFINEMENT
  void invQuantization(int* src, int* dst, int* si, int offsetX, int offsetY, int w, int h);
  void invQuantization(int* src, int* dst, int offsetX, int offsetY, int w, int h);
# else
  void invQuantization(int* src, int* dst, int* si, int w, int h);
  void invQuantization(int* src, int* dst, int w, int h);
# endif

private:
  // DCT
  void forward4x4(int** src, int** dst, int x, int y);

  // Inverse DCT
  void inverse4x4(int** src, int** dst, int x, int y);

  template <typename T>
  void idct4x4(int* src, T* dst, int x, int y, int w);

  // Quantization
  void quan4x4(int* src, int* dst, int x, int y, int w);

  // Inverse quantization
# if SI_REFINEMENT
  void invquan4x4(int* src, int* dst, int* si, int x, int y, int i, int j, int w);
  void invquan4x4(int* src, int* dst, int x, int y, int i, int j, int w);
# else
  void invquan4x4(int* src, int* dst, int* si, int x, int y, int w);
  void invquan4x4(int* src, int* dst, int x, int y, int w);
# endif

  //! Minimum mean square error
  /*! @param[out] x     De-quantized coefficient
      @param[in]  y     Side information
      @param[in]  z     Quantized coefficient
      @param[in]  step  Step size
      @param[in]  alpha Alpha value */
  void mmse(int* x, int* y, int* z, int step, double alpha);

private:
  const static int DctScaler[4][4];
  const static int DctShift[4][4];
  const static int IdctScaler[4][4];
  const static int IdctShift[4][4];

  Codec* _codec;
};

#endif // ENCODER_INC_TRANSFORM_H
