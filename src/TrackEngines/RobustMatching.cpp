#include "RobustMatching.h"

#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/gpu/gpu.hpp"
#include "opencv2/nonfree/gpu.hpp"
using namespace cv::gpu;

#include "opencv2/ocl/ocl.hpp"
#include "opencv2/nonfree/ocl.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"

using namespace cv;
using namespace cv::ocl;

std::vector<p_match> checkForLoop(cv::Mat &img_1_bw,cv::Mat &img_2_bw,Matrix3f &Homography)
{
	std::vector<cv::KeyPoint> keypoints1, keypoints2;
	std::vector<cv::DMatch> matches;	
	
	GpuMat img1, img2;
	img1.upload(img_1_bw);
	img2.upload(img_2_bw);
	
	float Hessian_threshold=800;
	SURF_GPU surf(Hessian_threshold);

	// detecting keypoints & computing descriptors
	GpuMat keypoints1GPU, keypoints2GPU;
	GpuMat descriptors1GPU, descriptors2GPU;
	surf(img1, GpuMat(), keypoints1GPU, descriptors1GPU);
	surf(img2, GpuMat(), keypoints2GPU, descriptors2GPU);

	// matching descriptors
	BFMatcher_GPU matcher(cv::NORM_L2);
	GpuMat trainIdx, distance;
	matcher.matchSingle(descriptors1GPU, descriptors2GPU, trainIdx, distance);

	// downloading results
	std::vector<float> descriptors1, descriptors2;
	surf.downloadKeypoints(keypoints1GPU, keypoints1);
	surf.downloadKeypoints(keypoints2GPU, keypoints2);
	surf.downloadDescriptors(descriptors1GPU, descriptors1);
	surf.downloadDescriptors(descriptors2GPU, descriptors2);
	
	BFMatcher_GPU::matchDownload(trainIdx, distance, matches);
	
	
	std::vector< cv::DMatch > good_matches;
	double max_dist = 0; double min_dist = 100;

	//-- Quick calculation of max and min distances between keypoints
	for( int i = 0; i < matches.size(); i++ )
	{ double dist = matches[i].distance;
	  if( dist < min_dist ) min_dist = dist;
	  if( dist > max_dist ) max_dist = dist;
	}


	for( int i = 0; i < matches.size(); i++ )
	{ if( matches[i].distance <= cv::max(2*min_dist, 0.02) )
	//{ if( matches[i].distance <= cv::max(2*min_dist, 0.2) )
	  { good_matches.push_back( matches[i]); }
	}


	
	std::vector<p_match> pmatches;
	
	for( size_t i = 0; i < good_matches.size(); i++ )
	{
	  
	    //-- Get the keypoints from the good matches
	    Point2f pt1= keypoints1[ good_matches[i].queryIdx ].pt ;
	    Point2f pt2= keypoints2[ good_matches[i].trainIdx ].pt ;
	    p_match newMatch(pt1.x,pt1.y,0,pt2.x,pt2.y,0);
	    pmatches.push_back(newMatch);
	}	
	
	if(pmatches.size()>4)
		Homography=HomographyFromMatchesRANSAC(pmatches,200);
	return pmatches;
}