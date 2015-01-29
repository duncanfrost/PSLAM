//Map tracking module.

#pragma once

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <sys/time.h>

#include <Eigen/Core>
using namespace Eigen;

#include "../AmoDefines.h"
#include "ImageProcess.h"
#include "../Primitives/obMap.h" 
#include "../Primitives/HomogeneousMatrix.h"
#include "GlobalTransfoEstim.h"
#include "../MapEngines/MapOptimiser.h" 
#include "../MapEngines/MapOptimiserEssential.h" 

//for fast corner detection
//#include <opencv/cvaux.h> //for opencv 2.4
#include <opencv2/opencv.hpp> // for later one


class MapTracker
{
public:
	MapTracker(){};
	~MapTracker(){};
	MapTracker(Camera *_cam,obMap *_map);
	void Init(Camera *_cam,obMap *_map);
	
	//estimate pose that minimizes reprojection error of map onto image
	//note that color image here is used only to define the color of the reconstructed point cloud for visualization
	void TrackFrame(cv::Mat &_current_img,cv::Mat *_current_img_col=NULL,bool useCol=false);
	
	//get current pose
	HomogeneousMatrix getPose(){return relPose*myMap->getKF(idrelKF)->getPose();};
	
	void updateRelPoseAfterMapOptim(){relPose=myMap->getKF(idrelKF)->getRelativePose();}
	
	//get id KF relative
	int getIdRelKF(){return idrelKF;};
	std::vector<int> &getClosestKFs(){return id_closestKF;};
private:
  
	//pointers o camera and map
	Camera *myCamera;
	obMap *myMap;
	
	//current relative pose between KFcurrent and KFclosest (in map)
	HomogeneousMatrix relPose;	
	
	int idrelKF;//id of KF relative pose is defined
	std::vector<int> id_closestKF;//id of current closest KF in map (active window)
	
	Matrix3f MotionPriorHomography;//can use MotionPriorHomography=Id or use constant mvt prior=> MotionPriorHomography=Homography(t-1)*Homography(t-2).inverse()
	//for instance in Kitti, when rotate, rotation is quiet smooth but large=> if we don not take motion prior into account lots of features are not matched
	
};

	
	
