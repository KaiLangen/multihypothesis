
#ifndef COMMON_INC_CONFIG_H
#define COMMON_INC_CONFIG_H

#define PI                      3.14159265
#define BLOCK_SIZE              4
#define MB_BLOCK_SIZE           16
#define DC_BITDEPTH             10

// Macros for encoder/decoder
#define WHOLEFLOW               1
#define AC_QSTEP                1
#define SKIP_MODE               1
#define MODE_DECISION           1
#define INTEGER_DCT             1
#define HARDWARE_FLOW           1
#define HARDWARE_LDPC           1
#define HARDWARE_CMS            1
#define HARDWARE_OPT            1
#define OBMC                    0

// Macros for encoder only
#ifdef ENCODER
# define DEBUG                  0
#endif

// Macros for decoder only
#ifdef DECODER
# define INVERSE_MATRIX         1
#endif

#include "types.h"

#endif // COMMON_INC_CONFIG_H

