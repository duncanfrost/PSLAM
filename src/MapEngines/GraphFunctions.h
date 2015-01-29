//few functions to work with graph structure of map

#pragma once


#include "../Primitives/obMap.h"
#define LevMarLikeConstant 1e-5;

struct kf_edgeNew
{
	int kf1,kf2;
	float weight;
};

struct scale_constraint
{
	int kf_id;
	float rescale;
};


	
//get edges in an innerWindow
std::vector<kf_edgeNew> getEdgesInInnerWin(obMap *myMap, std::vector<int> _innerWindowKFs);

