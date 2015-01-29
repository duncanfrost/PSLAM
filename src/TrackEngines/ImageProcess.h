#pragma once

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include "../AmoDefines.h"
#include "../CvWrapper.h"
#include "../Primitives/Camera.h"
#include "../MapEngines/GeoFunctions.h"

#include <Eigen/Core>
using namespace Eigen;

struct amoRect
{
 	amoRect(Vector2f _x1,Vector2f _x2){x1=_x1;x2=_x2;};
	
	Vector2f x1,x2;//top left and bottom right corners
	
	bool isInRect(Vector2f x){return (x[0]>=x1[0] && x[0]<x2[0] && x[1]>=x1[1] && x[1]<x2[1]);};
	bool isInRect(int x,int y){return (x>=x1[0] && x<x2[0] && y>=x1[0] && y<x2[1]);};
	bool isInRect(float x,float y){return (x>=x1[0] && x<x2[0] && y>=x1[0] && y<x2[1]);};
};

bool isInImage(Vector2i imageSize,Vector2f p);
bool isInImage(cv::Size imageSize,Vector2f p);
bool isInImage(cv::Mat &img,Vector2f p);
bool isInImageMargin(Vector2i imageSize,int margin,Vector2f p);
bool isInImageMargin(cv::Size imageSize,int margin,Vector2f p);
bool isInImageMargin(cv::Mat &img,int margin,Vector2f p);


//gradient of black and white image
void imageGradients(cv::Mat &Img,cv::Mat &gradxImg,cv::Mat &gradyImg);
//lower level image gradient computation functions
void GetGradX(cv::Mat &I, cv::Mat& I2);
void GetGradXnoOptim(cv::Mat &I, cv::Mat& dIx);//10x slower
inline double derivXF(BwImage & fr, int r, int c);
void GetGradY(cv::Mat &I, cv::Mat& dIy);
void GetGradYnoOptim(cv::Mat &I, cv::Mat& dIy);
double derivYF(BwImage & fr, int r, int c);


//get gaussian pyramid of image
void GetGauss(cv::Mat &I, cv::Mat& GI);
void GetGaussX(cv::Mat &I, cv::Mat& dIx);
double GaussXF(BwImage & fr, int r, int c);
void GetGaussY(cv::Mat &I, cv::Mat& dIy);
double GaussYF(BwImage & fr, int r, int c);

void extractPatch(cv::Mat &Img,Vector2f pixPos,cv::Mat &Patch,int patchSize=9);

float ScaleLevel(int _l);
float invScaleLevel(int _l);
//float ScaleCurrentToZeroLevel(int _l);
//float ScaleZeroToCurrentLevel(int _l);
double LevelNPos(double dRootPos, int nLevel);
double LevelZeroPos(double dLevelPos, int nLevel);
Vector2f ScaleCurrentToZeroLevel(Vector2f p,int _l);
Vector2f ScaleZeroToCurrentLevel(Vector2f p,int _l);

void drawLine(cv::Mat &imrgb,Vector2f p0,Vector2f p1,Vector3f color);
