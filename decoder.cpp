#include <cstdio>
#include <opencv2/opencv.hpp>

#include "types.h"
#include "obme.h"

using namespace cv;

int main()//int argc, char** argv)
{
  const char filename[] = "../bowing_cif.yuv";
  FILE* fp = fopen(filename, "rb");
  int height = 288;
  int width = 352;
  const int lsize = height * width;
  const int fsize = 3*(lsize>>1);
  int key_skip = 16;
  int nkeys = 6;
  int start = 80;
  int blocksize = 32;

  if (!fp)
    exit(EXIT_FAILURE);

  /* read all the key frames */
  fseek(fp, fsize * start, SEEK_SET);
  imgpel** refFrame = new imgpel * [nkeys];
  imgpel** refChromaU = new imgpel * [nkeys];
  imgpel** refChromaV = new imgpel * [nkeys];

  for(int i = 0; i < nkeys; i++) {
    fseek(fp, fsize * key_skip, SEEK_CUR);
    refFrame[i] = new imgpel[fsize];
    refChromaU[i] = new imgpel[lsize];
    refChromaV[i] = new imgpel[lsize];
    fread(refFrame[i], fsize, 1, fp);

    /* upsample the Chroma into new buffer */
    Mat mchromaU_in(height>>1, width>>1, CV_8UC1, refFrame[i] + lsize);
    Mat mchromaU_out(height, width, CV_8UC1, refChromaU[i]);
    resize(mchromaU_in, mchromaU_out, mchromaU_out.size(), 0, 0, INTER_CUBIC); 

    Mat mchromaV_in(height>>1, width>>1, CV_8UC1, refFrame[i] + 5*(lsize>>2));
    Mat mchromaV_out(height, width, CV_8UC1, refChromaV[i]);
    resize(mchromaV_in, mchromaV_out, mchromaV_out.size(), 0, 0, INTER_CUBIC); 
//    imshow("REAL Current Frame", mchroma2);
//    waitKey(0);
  }

  /* read the current frame */
  imgpel* currFrame = new imgpel[fsize];
  imgpel* currChromaU = new imgpel[lsize];
  imgpel* currChromaV = new imgpel[lsize];
  imgpel* currBgr = new imgpel[lsize * 3];
  mvinfo* mvs = new mvinfo[lsize / (blocksize * blocksize)];
  fseek(fp, fsize * (start + nkeys * key_skip / 2), SEEK_SET);
  fread(currFrame, fsize, 1, fp);
  
  /* get Chroma, upsample / blur */
  Mat mchromaU_in(height>>1, width>>1, CV_8UC1, currFrame + lsize);
  Mat mchromaU_out(height, width, CV_8UC1, currChromaU);
  Mat mchromaV_in(height>>1, width>>1, CV_8UC1, currFrame + 5*(lsize>>2));
  Mat mchromaV_out(height, width, CV_8UC1, currChromaV);
  resize(mchromaU_in, mchromaU_out, mchromaU_out.size(), 0, 0, INTER_CUBIC); 
  resize(mchromaV_in, mchromaV_out, mchromaV_out.size(), 0, 0, INTER_CUBIC); 
  imshow("REAL Current Chroma (V)", mchromaV_out);
  waitKey(0);
  imshow("REAL Current Chroma (U)", mchromaU_out);
  waitKey(0);

  Mat myuv(3*(height>>1), width, CV_8UC1, currFrame);
  Mat mbgr(height, width, CV_8UC3, currBgr);
  cvtColor(myuv, mbgr, COLOR_YUV2BGR_I420, 3);
  imshow("COLOUR Current Frame", mbgr);
  waitKey(0);

  /* Calculate MV's in U Chroma, then OBMC in the Luma domain */
  obme(refChromaU, currChromaU,
       refChromaV, currChromaV, mvs,
       nkeys, width, height, blocksize);

  memset(currFrame, 0x0, lsize);
  obmc(refFrame, currFrame, mvs, width, height, blocksize);
  Mat myuv2(3*height/2, width, CV_8UC1, currFrame);
  Mat mbgr2(height, width, CV_8UC3, currBgr);
  cvtColor(myuv2, mbgr2, COLOR_YUV2BGR_I420, 3);
  imshow("ESTIMATED Current Frame", mbgr2);
  waitKey(0);

  delete [] currFrame;
  delete [] currChromaU;
  delete [] currChromaV;
  delete [] currBgr;
  delete [] mvs;
  for(int i = 0; i < nkeys; i++) {
    delete [] refFrame[i];
    delete [] refChromaU[i];
    delete [] refChromaV[i];
  }
  delete [] refFrame;
  delete [] refChromaU;
  delete [] refChromaV;
  fclose(fp);

  return 0;
}
