//sample to show how robust matching is performed using ORB
//press space bar to define new reference image

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"



#include <Eigen/Core>
#include <vector_types.h>
using namespace Eigen;

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);
void addDrawFunction(void) ;
//std::vector<Vector2i> GetListCorner(cv::Mat &Img);

//Visualization
VisualisationModule *VisuEngine;
//image acquisition
VideoSourceLiveCV *myVideoSource;
Camera myCamera;//camera object, calibration ...

cv::Mat ImgDisplay;
cv::Mat ImgCurrent;
cv::Mat ImgFirst;

//ORB things:
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

const int LOOP_NUM = 1;//was there to compute timing
const int GOOD_PTS_MAX = 50;
const float GOOD_PORTION = 0.15f;


int main(int argc, char** argv)
{
	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->setLvlAcquisition(1);
	myVideoSource->initCam();
	myCamera=*myVideoSource->getPointerCamera();

	myVideoSource->grabNewFrame();
	cv::Mat ImCam = myVideoSource->GetFramePointer()->clone();
	ImgDisplay.create(ImCam.size().height, ImCam.size().width*2, CV_8UC3);
	
	//init display
	uchar3 black;black.x=0;black.y=0;black.y=0;
	for(int y=0;y<ImgDisplay.size().height;y++)
		for(int x=0;x<ImgDisplay.size().width;x++)
			ImgDisplay.at<uchar3>(y, x)=black;

	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("SurfMatching",&ImgDisplay);
	VisuEngine->setOnKeyPress(&processNormalKeys);
	VisuEngine->setOnDraw(&addDrawFunction);

	cv::gpu::printShortCudaDeviceInfo(cv::gpu::getDevice());
		
	std::cout<<"press space to select first image"<<std::endl;

	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


struct matchDisp
{
	Vector2f pt1;
	Vector2f pt2;
};

std::vector<matchDisp> listMatchesDisp;
std::vector<Vector2f> scene_corners_disp(4);

void Idle(void) 
{
	myVideoSource->grabNewFrame();
	ImgCurrent=myVideoSource->GetFramePointer()->clone();

	amoTimer timerLoop;timerLoop.start();

	if(!ImgFirst.empty())//there is a ref to do matching
	{
		
		std::vector<cv::KeyPoint> keypoints1, keypoints2;
		std::vector<cv::DMatch> matches;	
		
		cv::Mat &cpu_img2=ImgCurrent;
		cv::Mat &cpu_img1=ImgFirst;

		//matchImagesWithSurfOcl(cpu_img1,cpu_img2,keypoints1,keypoints2,matches);
		
		
		cv::Mat img_2_bw;cv::cvtColor(ImgCurrent,img_2_bw,CV_RGB2GRAY);
		cv::Mat img_1_bw;cv::cvtColor(ImgFirst,img_1_bw,CV_RGB2GRAY);
		
		amoTimer timer;timer.start();
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

		//std::cout << "FOUND " << keypoints1GPU.cols << " keypoints on first image" << std::endl;
		//std::cout << "FOUND " << keypoints2GPU.cols << " keypoints on second image" << std::endl;

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
		timer.stop("Extraction");
		amoTimer timer2;timer2.start();
		
		BFMatcher_GPU::matchDownload(trainIdx, distance, matches);
		
		timer2.stop("matching");

		
		std::vector< cv::DMatch > good_matches;
		if(1)//based on distance
		{
			double max_dist = 0; double min_dist = 100;

			//-- Quick calculation of max and min distances between keypoints
			for( int i = 0; i < matches.size(); i++ )
			{ double dist = matches[i].distance;
			  if( dist < min_dist ) min_dist = dist;
			  if( dist > max_dist ) max_dist = dist;
			}
			//-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist,
			//-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
			//-- small)
			//-- PS.- radiusMatch can also be used here.

			for( int i = 0; i < matches.size(); i++ )
			{ if( matches[i].distance <= cv::max(2*min_dist, 0.02) )
			//{ if( matches[i].distance <= cv::max(2*min_dist, 0.2) )
			  { good_matches.push_back( matches[i]); }
			}
		}
		else //keep only best GOOD_PTS_MAX
		{
			//-- Sort matches and preserve top 10% matches
			std::sort(matches.begin(), matches.end());
			double minDist = matches.front().distance,
			      maxDist = matches.back().distance;

			const int ptsPairs = std::min(GOOD_PTS_MAX, (int)(matches.size() * GOOD_PORTION));
			for( int i = 0; i < ptsPairs; i++ )
			{
			    good_matches.push_back( matches[i] );
			}
		}
		
		 //display ref and matches
		  for(int y=0;y<ImgCurrent.size().height;y++)
			  for(int x=0;x<ImgCurrent.size().width;x++)
				  ImgDisplay.at<uchar3>(y, x)=ImgFirst.at<uchar3>(y, x);
			  
		//create matches
		listMatchesDisp.clear();
		for(int m=0;m<good_matches.size();m++)
		{
			matchDisp newMatchDisp;
			int i1 = good_matches[m].queryIdx;
			int i2 = good_matches[m].trainIdx;
			assert(i1 >= 0 && i1 < static_cast<int>(keypoints1.size()));
			assert(i2 >= 0 && i2 < static_cast<int>(keypoints2.size()));

			const cv::KeyPoint &kp1 = keypoints1[i1], &kp2 = keypoints2[i2];
			newMatchDisp.pt1=Vector2f(kp1.pt.x,kp1.pt.y);
			newMatchDisp.pt2=Vector2f(kp2.pt.x,kp2.pt.y);
			
			listMatchesDisp.push_back(newMatchDisp);
		}


		//-- Localize the new image
		std::vector<Point2f> obj;
		std::vector<Point2f> scene;

		for( size_t i = 0; i < good_matches.size(); i++ )
		{
		    //-- Get the keypoints from the good matches
		    obj.push_back( keypoints1[ good_matches[i].queryIdx ].pt );
		    scene.push_back( keypoints2[ good_matches[i].trainIdx ].pt );
		}		
		//get homography
		//-- Get the corners from the image_1 ( the object to be "detected" )
		std::vector<Point2f> obj_corners(4);
		obj_corners[0] = cvPoint(0,0);
		obj_corners[1] = cvPoint( cpu_img1.cols, 0 );
		obj_corners[2] = cvPoint( cpu_img1.cols, cpu_img1.rows );
		obj_corners[3] = cvPoint( 0, cpu_img1.rows );
		
		std::vector<Point2f> scene_corners(4);
		if(obj.size()>=4)
		{
			Mat H = findHomography( obj, scene, CV_RANSAC );
			perspectiveTransform( obj_corners, scene_corners, H);
		}
		
		for(int i=0;i<4;i++)
			scene_corners_disp[i]=Vector2f(scene_corners[i].x,scene_corners[i].y);
		
		std::cout<<obj.size() <<"matches"<<std::endl;
		 
	}

	//display input image in right of window
	for(int y=0;y<ImgCurrent.size().height;y++)
		for(int x=0;x<ImgCurrent.size().width;x++)
			ImgDisplay.at<uchar3>(y, x+ImgCurrent.size().width)=ImgCurrent.at<uchar3>(y, x);

	timerLoop.stop();
	
	//timer.stop();
	VisuEngine->drawWindows();
}
void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			break;
		case ' '://esc
			ImgFirst=ImgCurrent.clone();
			break;
	}


	glutPostRedisplay();
}


void addDrawFunction(void) 
{	
	
	//get max response feature
	
	set2DGLProjection();
	glPointSize(3.0);
	//glColor3f(0,1,0);
	cv::RNG& rng=cv::theRNG();
	glLineWidth(2.);
	for(int i=0;i<listMatchesDisp.size();i++)
	{
		glColor3f(rng(256)/256.,rng(256)/256.,rng(256)/256.);
		glBegin(GL_LINES);
		glVertex2f(listMatchesDisp[i].pt1[0],listMatchesDisp[i].pt1[1]);
		glVertex2f(listMatchesDisp[i].pt2[0]+ImgFirst.size().width,listMatchesDisp[i].pt2[1]);
		glEnd();
	}
	
	glLineWidth(4.);
	glColor3f(0,1.,0);
	for(int i=0;i<4;i++)
	{
		glBegin(GL_LINES);
		glVertex2f(scene_corners_disp[i][0]+ImgFirst.size().width,scene_corners_disp[i][1]);
		glVertex2f(scene_corners_disp[(i+1)%4][0]+ImgFirst.size().width,scene_corners_disp[(i+1)%4][1]);
		glEnd();
	}
	
	glColor3f(1,1,1);
	unset2DGLProjection();
}


/*
#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/features2d.hpp"

using namespace cv;

void readme();


int main( int argc, char** argv )
{
  if( argc != 3 )
  { readme(); return -1; }

  Mat img_1 = imread( argv[1], CV_LOAD_IMAGE_GRAYSCALE );
  Mat img_2 = imread( argv[2], CV_LOAD_IMAGE_GRAYSCALE );

  if( !img_1.data || !img_2.data )
  { std::cout<< " --(!) Error reading images " << std::endl; return -1; }

  //-- Step 1: Detect the keypoints using SURF Detector
  int minHessian = 400;

  SurfFeatureDetector detector( minHessian );

  std::vector<KeyPoint> keypoints_1, keypoints_2;

  detector.detect( img_1, keypoints_1 );
  detector.detect( img_2, keypoints_2 );

  //-- Step 2: Calculate descriptors (feature vectors)
  SurfDescriptorExtractor extractor;

  Mat descriptors_1, descriptors_2;

  extractor.compute( img_1, keypoints_1, descriptors_1 );
  extractor.compute( img_2, keypoints_2, descriptors_2 );

  //-- Step 3: Matching descriptor vectors using FLANN matcher
  FlannBasedMatcher matcher;
  std::vector< DMatch > matches;
  matcher.match( descriptors_1, descriptors_2, matches );

  double max_dist = 0; double min_dist = 100;

  //-- Quick calculation of max and min distances between keypoints
  for( int i = 0; i < descriptors_1.rows; i++ )
  { double dist = matches[i].distance;
    if( dist < min_dist ) min_dist = dist;
    if( dist > max_dist ) max_dist = dist;
  }

  printf("-- Max dist : %f \n", max_dist );
  printf("-- Min dist : %f \n", min_dist );

  //-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist,
  //-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
  //-- small)
  //-- PS.- radiusMatch can also be used here.
  std::vector< DMatch > good_matches;

  for( int i = 0; i < descriptors_1.rows; i++ )
//  { if( matches[i].distance <= max(2*min_dist, 0.02) )
  { if( matches[i].distance <= max(2*min_dist, 0.2) )
    { good_matches.push_back( matches[i]); }
  }

  //-- Draw only "good" matches
  Mat img_matches;
  drawMatches( img_1, keypoints_1, img_2, keypoints_2,
               good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
               vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );

  //-- Show detected matches
  imshow( "Good Matches", img_matches );

  for( int i = 0; i < (int)good_matches.size(); i++ )
  { printf( "-- Good Match [%d] Keypoint 1: %d  -- Keypoint 2: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx ); }

  waitKey(0);

  return 0;
}


void readme()
{ std::cout << " Usage: ./SURF_FlannMatcher <img1> <img2>" << std::endl; }*/
