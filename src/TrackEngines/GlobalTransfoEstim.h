#pragma once

//stuff to warp images, find relative transformation between them etc...

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include "../AmoDefines.h"
#include "../CvWrapper.h"
#include "../Primitives/Camera.h" 
#include "../MapEngines/GeoFunctions.h"
#include "../TrackEngines/ImageProcess.h" 
#include "../Primitives/HomogeneousMatrix.h" 

#include <Eigen/Core>
using namespace Eigen;

#include "omp.h"

cv::Mat warpImageInv(cv::Mat &Img,Matrix3f _homography);
//check how much of two images overlap when they are linked with _homography, modulo=sample each modulo pixel to compute overlap
float getOverlapFromHomography(cv::Size imageSize,Matrix3f _homography,short modulo);

//apply warps to image points in pixels, result in pixel too
Vector2f applyRotation(Vector2f pos,Vector3f w,Camera* _cam);
Vector2f applyRotationPyr(Vector2f pos,Vector3f w,Camera* _cam,int _lvl);//apply same transformation but in lvl of pyramid
Matrix2f getSO2RotationFromAplha(float _a);
Vector3f invSE2(Vector3f p);
Vector2f applySE2(Vector2f pos,Vector3f p,Camera* _cam);
Vector2f applySE2Pyr(Vector2f pos,Vector3f p,Camera* _cam,int _lvl);
Vector3f SE2toSO3(Vector3f _p,Camera* _cam);

//transform images
cv::Mat translateImage(cv::Mat &Img,Vector2f t);
cv::Mat rotateImage(cv::Mat &Img,Vector3f w,Camera* _cam);
cv::Mat applySE2Image(cv::Mat &Img,Vector3f p,Camera* _cam);

//find transformation between two images

//estimate translation from T to I; so that T(x)= I(w(x))
Vector2f estimateTranslation(cv::Mat &T,cv::Mat &I,int modulo=1);

//estimate SE2 from T to I; so that T(x)= I(w(x))
Vector3f estimateSE2(cv::Mat &T,cv::Mat &I, Camera* _cam,int modulo=1);
Vector3f estimateSE2Pyr(cv::Mat &T,cv::Mat &I, Camera* _cam,int pyr_lvl=0,int modulo=1,int _max_iter=50);

//estimate rotation from T to I; so that T(x)= I(w(x))
Vector3f estimateRotation(cv::Mat &T,cv::Mat &I, Camera* _cam,int modulo=1);
Vector3f estimateRotationInverse(cv::Mat &T,cv::Mat &I, Camera* _cam,int modulo=1);//same but inverse appproach= gradients computed on T
Vector3f estimateRotationPyr(cv::Mat &T,cv::Mat &I, Camera* _cam,int pyr_lvl=0,int modulo=1);


//not really tz here but more scale s used as (1+s)Rx
Vector2f applySE2Tz(Vector2f pos,Vector4f w,Camera* _cam);
Vector2f applySE2TzPyr(Vector2f pos,Vector4f w,Camera* _cam,int _lvl);//apply same transformation but in lvl of pyramid
cv::Mat applySE2TzImage(cv::Mat &Img,Vector4f p,Camera* _cam);
HomogeneousMatrix SE2TztoSO3(Vector4f _p,float mean_depth,Camera* _cam);//SE3 + translation along z (for relocalization)
//SE3 + translation along z (for relocalization) return normalized final_residue
float estimateTzPyr(cv::Mat &T,cv::Mat &I, Camera* _cam,int pyr_lvl=0,int modulo=1,int _max_iter=50);
float estimateSE2TzPyr(cv::Mat &T,cv::Mat &I, Camera* _cam,Vector4f &res,int pyr_lvl=0,int modulo=1,int _max_iter=50);

#include "../VisuOdo/motionEstim.h"

Vector3f estimateRotation(std::vector<p_match> &matches, Camera* _cam);
Vector3f estimateRotationInv(std::vector<p_match> &matches, Camera* _cam);

//rotate image but resize output so that all image is still in there=> create an offset
cv::Mat rotateImageKeepAll(cv::Mat &Img,Vector3f w,Camera* _cam,Vector2f &offsetRes);
