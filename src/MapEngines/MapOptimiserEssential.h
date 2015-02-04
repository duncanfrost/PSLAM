//bundle adjustment like functions: 
//optimise xEx for up to scale transformation
//and optimise over the variance of the local 3d features for scale

#pragma once


#include "../Primitives/obMap.h"
#define LevMarLikeConstant 1e-5;


class MapOptimiserEssential
{
public:
	MapOptimiserEssential(obMap *_myMap);
	~MapOptimiserEssential(){};
	
	//will optimise scale and relative positions of KFs
	//if innerWindow corresponds to all the frames of the submap, ie no
	//link to fix KF then will have to choose one KF of innerwindow as fix
	void optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	void optimiseInnerWindowRobust(std::vector<int> &_innerWindowKFs,int nb_iter=10);

private:
	//map to be optimised
	obMap *myMap;
	float gain;
	
	//KFs in inner window => with scale and relative poses to be optimised
	std::vector<int> innerWindowKFs;
};
