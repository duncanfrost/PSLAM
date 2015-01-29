
#pragma once

#include "motionEstim.h"

#include "../AmoDefines.h"
#include "../TrackEngines/ImageProcess.h"
#include "../MapEngines/HomogFinder.h"
#include "../TrackEngines/GlobalTransfoEstim.h"
using namespace Eigen;

#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

class matcherWrapper{

public:

	matcherWrapper();	
	// deconstructor
	~matcherWrapper(){};
	
	//main tracking functions
	void InitWithRef(cv::Mat *_img_p);
	void match(cv::Mat *_img_p);
	void useKLT(bool _b){useKLTpoints=_b;};
	void useKLTonly(bool _b){useKLTpointsOnly=_b;};
  
	Matrix3f getHomography(){return Homography;};
	void setHomography(Matrix3f _h){Homography=_h;};
	
	std::vector<p_feat> getFeaturesp(){return matcher.getFeaturesp();};
	std::vector<p_match> getCurrentMatches(){return matches;};
	std::vector<p_match> getCurrentMatchesBucket(){return matchesBucket;};

private:
	int mlvl_viso;
	Matcher matcher;
	
	//klt stuff
	bool useKLTpoints;
	bool useKLTpointsOnly;
	int MAX_COUNT_KLT;
	int dist_klt_feat;
	
	
	//output:
	Matrix3f Homography;//homog at lvl 0
	std::vector<p_match> matches;
	std::vector<p_match> matchesBucket;


};


