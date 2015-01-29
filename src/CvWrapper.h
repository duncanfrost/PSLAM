//definition of runnable object
#ifndef __CVWRAP_H___
#define __CVWRAP_H___

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "AmoDefines.h"

#include <Eigen/Core>
using namespace Eigen;

template<class T> class Image
{
  private:
  IplImage* imgp;
  IplImage imgp0;
  
  public:
  Image(IplImage* img=0) {imgp=img;width=img->width;}
  Image(cv::Mat &img) {imgp0=img;imgp=&imgp0;width=img.size().width;}

  ~Image(){imgp=0;}
  void operator=(IplImage* img) {imgp=img;}
  inline T* operator[](const int rowIndx) {return ((T *)(imgp->imageData + rowIndx*imgp->widthStep));}
  
    int width;
};


typedef struct{
  unsigned char b,g,r;
} RgbPixel;

typedef struct{
  float b,g,r;
} RgbPixelFloat;

typedef Image<RgbPixel>       RgbImage;
typedef Image<RgbPixelFloat>  RgbImageFloat;
typedef Image<unsigned char>  BwImage;
typedef Image<float>          BwImageFloat;

float getColorSubpix(const cv::Mat& img, Vector2f pt);
float getColorSubpixf(const cv::Mat& img, Vector2f pt);
cv::Vec3b getColorSubpixRGB(const cv::Mat& img, cv::Point2f pt);

cv::Mat getDisplayableImage(cv::Mat& img);

#endif
