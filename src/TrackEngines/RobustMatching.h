#pragma once

//list of functions to help for some robust detection
//not used yet.

#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <objectr3d/AmoDefines.h>
#include <objectr3d/CvWrapper.h>
#include <objectr3d/CvWrapper.h>
#include <objectr3d/Camera.h>
#include <objectr3d/GeoFunctions.h>
#include <objectr3d/motionEstim.h>
#include <objectr3d/HomogFinder.h>

#include <Eigen/Core>
using namespace Eigen;

//extract ORB features for images, match them and return matches and homography between them
//note that it has been coded to test robust matching but if used later, we ll have to 
//pre extract ORB features in each keyframe so that only the matching will have to be done to check for possible loop closures
std::vector<p_match> checkForLoop(cv::Mat &img1,cv::Mat &img2,Matrix3f &Homography);