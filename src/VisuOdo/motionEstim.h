
#ifndef MOTION_ESTIM_H
#define MOTION_ESTIM_H

#include "viso.h"



 
  
class MotionEstimation{

public:

	// constructor, takes as inpute a parameter structure
	MotionEstimation (float fu,float fv,float cu,float cv);
	
	// deconstructor
	~MotionEstimation ();
  
	bool estimateMotion (std::vector<p_match> p_matched);
	bool estimateMotionAndRemoveOutliers (std::vector<p_match> &p_matched);
	void getEstimatedTransfo(float *rot,float *trans);


	private:
	Viso::Matrix               smallerThanMedian (Viso::Matrix &X,double &median);
	bool                 normalizeFeaturePoints (std::vector<p_match> &p_matched,Viso::Matrix &Tp,Viso::Matrix &Tc);
	void                 fundamentalMatrix (const std::vector<p_match> &p_matched,const std::vector<int32_t> &active,Viso::Matrix &F);
	void                 EtoRt(Viso::Matrix &E,Viso::Matrix &K,std::vector<p_match> &p_matched,Viso::Matrix &X,Viso::Matrix &R,Viso::Matrix &t);
	int32_t              triangulateChieral (std::vector<p_match> &p_matched,Viso::Matrix &K,Viso::Matrix &R,Viso::Matrix &t,Viso::Matrix &X);
	void                 EtoRtOutlierReject(Viso::Matrix &E,Viso::Matrix &K,std::vector<p_match> &p_matched,Viso::Matrix &X,Viso::Matrix &R,Viso::Matrix &t);
	int32_t              triangulateChieralOutlierReject (std::vector<p_match> &p_matched,Viso::Matrix &K,Viso::Matrix &R,Viso::Matrix &t,Viso::Matrix &X,Viso::Matrix &isXoulier);
	std::vector<int32_t> getInlier (std::vector<p_match> &p_matched,Viso::Matrix &F);
	std::vector<int32_t> getRandomSample(int32_t N,int32_t num);

	//cam params
	float fu,fv,cu,cv;

	// parameters
	float pitch          ;
	float ransac_iters   ;
	std::vector<int32_t>           inliers;    // inlier set
	float inlier_threshold;
	float motion_threshold;

	//resulting motion
	Viso::Matrix R,t;

};

#endif // MOTION_ESTIM_H

