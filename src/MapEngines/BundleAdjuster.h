#pragma once

#include "../Primitives/Camera.h"
#include "../TrackEngines/ImageProcess.h"
#include "../Primitives/HomogeneousMatrix.h"
#include "../Primitives/obMap.h"

#define LevMarLikeConstant 1e-5;


struct BAcamera
{
	HomogeneousMatrix pose;//current pose
	HomogeneousMatrix poseNew;//updated pose, copied to pose if improve reproj error
	
	int nb_measures;//number of points of bundle visible in cam
	
	bool fixed;//is the camera position to be updated //is used to force cam to be fixed
	
	int id_optim;//position in update vector
	int id_main;//position in keyframe list in main map
};

struct BApoint
{
	Vector3f position;//current pose
	Vector3f positionNew;//updated pose, copied to pose if improve reproj error
	float confidence;
	
	int nb_measures;//number of cameras the point is measured in

	int id_optim;//position in update vector
	int id_main;//position in point list in main map
	int id_kf_main;//position in point list in main map
};

struct BAmeasure
{
	short kf_id;//id of corresponding kf in mCams
	short pt_id;//id of corresponding point in mPoints
	Vector2f coord;//2d meter position in lvl 0 pyr
	char lvl_meas;
	
	short id_feat_main;//position of measure in feature list of kf in map //only used for outlier rejection
};

struct OutlierMeasure
{
	short kf_id_main;//id of corresponding kf in mCams
	char lvl_meas;
	short id_feat_main;//position of measure in feature list of kf in map
};

enum ConvergenceResult { Converged, NotConverged, Diverged };

/*!
 * This implements the local bundle adjustment that happens in keyframes
 */

class BundleAdjuster
{
public:
	BundleAdjuster(obMap *_myMap=0);
	void optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter=10,bool robust=true);
	
	
    ~BundleAdjuster(){}
#ifdef USE_OMP_C	
	void setMoreImportantStuffToDo(bool *_moreImportantStuffWaiting,omp_lock_t *_lock_check_more_prior);
#endif	
	//fill Bundle
	void addCamera(HomogeneousMatrix _pose,int _id,bool _fixed);
	void addPoint(Vector3f _position,int _id_kf_pt,int _id,float _recAngle);
	//addMeasure should be called after addCamera and addPoint
	void addMeasure(short _idKFmeas,short _idKFPoint,short _idPoint,Vector2f coord,char _lvl,short id_feat);
	
	//run iteration of BA or point optimisation
	void optimise(int nb_iter=10,bool robust=true);
	//run iteration of BA or point optimisation
	//one Tukey coef for all measures
	void BundleAdjustRobust(int nb_iter=10);
	//one Tukey coef per Keyframe
	void BundleAdjustRobust2(int nb_iter=10);
	//run iteration of BA or point optimisation
	void OptimisePointPosition(int nb_iter=10);
	
	//outlier management
	//put measure in oulier list, if associated point does not have enough measure anymore=> remove from to optim list
	void rejectMeasure(int m);
	//put point in oulier list => remove all associated measures and keep them as outliers
	void rejectPoint(int p);
	//do not optimse point any more, measures are kept for cam updates only
	void doNotOptimisePoint(int p);
	
	//id management, get Point id in mPoints from its idOptim
	int IdOptimToIdBA(int i);
	
	//return updated value
	HomogeneousMatrix getUpdatedCamPose(int id_main);
	Vector3f getUpdatedPointPosition(int id_kf_main,int id_main);
	
	//return resulting outliers
    int nbOultiers(){return mOutliers.size();}
    OutlierMeasure &getOultiers(int i){return mOutliers[i];}
	
	//make it clear functions :
	void BundleAdjust(int nb_iter=10,bool LM=false);//like BundleAdjustRobust without Tukey function
	void BundleAdjustNoOptim(int nb_iter);//like BundleAdjustRobust without optimisation proposed by [Engels in bundle adjustment rules]

	bool hasBAConverged(){return hasConverged;}	
	bool hasBAGoneWrong(){return goneWrong;}
    ConvergenceResult getConvergenceResult();
private:
	//for io with Map
	obMap *myMap;
	float gain;
  
  
	//for BA
	Camera *myCamera;
	
	//list of cameras
	std::vector<BAcamera> mCams;
	int nbCamsToUpdate;
	
	//list of points
	std::vector<BApoint> mPoints;
	int nbPointsToUpdate;
		
	//list of measures
	std::vector<BAmeasure> mMeasures;
	
	//list of resulting outliers
	std::vector<OutlierMeasure> mOutliers;
	
	//convergence test
	bool hasConverged;
	//threshold of point position update
	float point_translation_small;
	//threshold of cam position translation update
	float cam_translation_small;
	
	//check for weird behavior
	bool goneWrong;
	
	//interruption variables (if more important task to be done)
#ifdef USE_OMP_C
	bool canBeInterrupted;
	omp_lock_t *lock_check_more_prior;
	bool *moreImportantStuffWaiting;
#endif
	
};



