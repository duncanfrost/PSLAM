//bundle adjustment like functions: optimise over the variance of the 
//local 3d features

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


/*!
 * \brief The MapOptimiser class implements the global map point minimisation.
 */
class MapOptimiser
{
public:
	MapOptimiser(obMap *_myMap);
    ~MapOptimiser(){}
	
	//use only keyframe to keyframe matches(no mapPoints=> good for marginalisation)
	void optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	//use mapPoints structure: minimise variance of point with all its corresponding local features
	void optimiseInnerWindow2(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	void optimiseInnerWindowRobust(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	void optimiseInnerWindowRobust2(std::vector<int> &_innerWindowKFs,int nb_iter=10);
	
	//do optimisation on scale: want scales on each side of an edge to be constant (depending on score edge)
	//and want scale to respect constaints
	void optimiseScale(std::vector<int> &_innerWindowKFs,std::vector<scale_constraint> &mScaleConstraints,int nb_iter=50);
	
	//readjust map point positions
	void optimiseMapPoints(std::vector<int> &_innerWindowKFs);
	
	//get optimal relative pose and scale between two keyframes and information matrix resulting for the alignment of theor local features
	void getRelativePoseAndScale(int _kf1,int _kf2,float &optRelScale,HomogeneousMatrix &optRelPose,float &infoScale,MatrixXf &infoPose);

    std::vector<kf_edgeNew> getEdgesInInnerWin(obMap *myMap, std::vector<int> _innerWindowKFs);


private:
	//map to be optimised
	obMap *myMap;
	float gain;
	
	//KFs in inner window => with scale and relative poses to be optimised
	std::vector<int> innerWindowKFs;
};
