//Detection and posterior segmentation
//Object Detection defines the detection primitives using a
//part based detector, create color models and estimate pixel wise foreground background proba and 
//defines the output structure

//TODO : needs a bit of cleaning

#pragma once

#include <objectr3d/Camera.h>
#include <objectr3d/MapPoint.h>
#include <objectr3d/KeyFrame.h>
#include <objectr3d/obMap.h>

#define LevMarLikeConstant 1e-5;

struct PGOcamera
{
	HomogeneousMatrix pose;//current pose
	
	int id_optim;//position in update vector, -1 if fixed
	int id_main;//position in keyframe list in main map
	float scale;
};
struct PGOedge
{
	int id_cam1;//position in mCams of first cam
	int id_cam2;//position in mCams of second cam
	
	HomogeneousMatrix pose12;//current relative pose from coord in 1 to coord in 2
	MatrixXf InfoMatrix;//current relative pose
	//for sim3
	float scale12;//scale to be applied to 1 to match 2
	float InfoScale;
};

class PoseGraphOptimiser
{
public:
	PoseGraphOptimiser(obMap *_myMap);
	void optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	~PoseGraphOptimiser(){};
	
	//fill Bundle
	void addCamera(HomogeneousMatrix _pose,int _id,bool _fixed);
	void addEdge(int _id1,int _id2,HomogeneousMatrix pose12);
	void addEdge(int _id1,int _id2,HomogeneousMatrix pose12,MatrixXf _InfoMatrix);
	void addEdge(int _id1,int _id2,float scale12,float _infoScale,HomogeneousMatrix pose12,MatrixXf _InfoMatrix);
	
	//run iteration of BA or point optimisation
	void optimise(int nb_iter=10);
	void optimiseSim3(int nb_iter=10);
	
	//return updated value
	HomogeneousMatrix getUpdatedCamPose(int id_main);
	float getUpdatedCamScale(int id_main);

	bool hasBAConverged(){return hasConverged;}	
private:
	//for io with Map
	obMap *myMap;
	float gain;
  
  
	//for PGO
	//list of cameras
	std::vector<PGOcamera> mCams;
	int nbCamsToUpdate;
	
	//list of edges
	std::vector<PGOedge> mEdges;
	int nbPointsToUpdate;
	
	//convergence test
	bool hasConverged;
	float cam_translation_small;

};




