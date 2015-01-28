//Detection and posterior segmentation
//Object Detection defines the detection primitives using a
//part based detector, create color models and estimate pixel wise foreground background proba and 
//defines the output structure

//TODO : needs a bit of cleaning

#pragma once

#include <algorithm>

#include <objectr3d/AmoDefines.h>
#include <objectr3d/KeyFrame.h>
#include <objectr3d/Camera.h>
//#include <objectr3d/BundleAdjuster.h>
//#include <objectr3d/PoseGraphOpt.h>

#include "opencv2/highgui/highgui.hpp"

#define KFoverlap 0.8//maximum overlap between two frames if rec angle is small
#define KFoverlapMin 0.5//minimum overlap between two KF
#define KFreconstructionAngle 10*3.14/180//10 deg //average reconstruction angle of points between two keyframe considered good to create new KF
//#define KFreconstructionAngle 2*3.14/180//10 deg //average reconstruction angle of points between two keyframe considered good to create new KF


class obMap
{
public:
	obMap();
	~obMap(){};
	
	//get camera pointer
	void InitCam(Camera *_myCamera){myCamera=_myCamera;};
	Camera *getCamera(){return myCamera;};
	
	//KF control
	int getNbKeyFrames(){return KeyFrameList.size();};
	KeyFrame* getKF(int i){return &KeyFrameList[i];}
	
	//if want to start a new submap (ie part of map with no known 
	//relative pose and scale with existing one); return id of new KF create
	int createUnconnectedNewKF(cv::Mat &_img);
	
	//if map tracker need a new key frame=> provide current image,relative pose with one KF and current neigbors
	//create new key frame and update neigbors with new KF
	void createNewKeyFrame(cv::Mat &_current_img,HomogeneousMatrix &relPose,int &idrelKF,std::vector<int> &id_closestKF);
	//for simulation
	void createNewKeyFrame(cv::Mat &_current_img,HomogeneousMatrix &_pose);
	void createNewEdge(int kfIdFrom,int kfIdTo);
	
	
	//function to access graph of connected Keyframes
	//get the direct neigbors of a set of keyframes 
	std::vector<int> getDirectNeigbors(std::vector<int> &_KFs);
	//get the direct neigbors of a set of keyframes and check using estimated homography if should be added to neigbor of current frame
	std::vector<int> getDirectNeigborsWithGoodEstimatedOverlap(std::vector<int> &_KFs);
	std::vector<int> getDirectNeigborsWithGoodEstimatedOverlapDepth(std::vector<int> &_KFs,int _d);
	//get submap connected with one KF, ie list of KF with distance < depth (result includes input)
	std::vector<int> getConnectedKeyframes(int _idkf,int _depth=-1){std::vector<int> _idkfv;_idkfv.push_back(_idkf);return getConnectedKeyframes(_idkfv,_depth);};
	std::vector<int> getConnectedKeyframes(std::vector<int> _idkf,int _depth=-1);
	
	//get information matrix corresponding to min variance optimisation between local features of two neigboring keyframes
	//can be used in PGO later
	MatrixXf getInformationMatrixMatches(int kfOrig,NeigbourKFNew &neigb, Camera *myCamera);
	void getRelativePoseAndScale(int _kf1,int _kf2,float &optRelScale,HomogeneousMatrix &optRelPose,float &infoScale,MatrixXf &infoPose,int nb_iter=20);

	//if a KF in map had a better fundamental matrix => new local features must have been create
	//=>new map points can be created
	void checkNeigborsForPoints(int idKF,bool recheckMatchedFeatures=false);
	
	//io functions
	void saveToStream(std::ofstream &fout);
	void loadFromStream(std::ifstream &fout);	
	void saveToFile(const char *filename);
	void loadFromFile(const char *filename);
	
	//remove unused points from KF (usually they are the results of duplicated points)
	void removeUnusedPoints(int _kfId);
	
	int getNbMapPoints();
	int getNbUsedMapPoints();
	
	//check if Kf idKF is matching with one of its neigbors distant of depth
	void checkForSmallLoop(int _idKF,int _depth_min,int _depth_max);
	
	//get reprojection error
	float getReprojectionError();
	

protected:
	//camera pointer
	Camera *myCamera;
	
	//map
	std::vector<KeyFrame> KeyFrameList;

};
