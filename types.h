#ifndef TYPES_H_
#define TYPES_H_

typedef char imgpel;

struct mvinfo {
  int iCx;
  int iCy;
  int iMvx;
  int iMvy;
  unsigned int SAD; 
  int frame;
};

#endif
