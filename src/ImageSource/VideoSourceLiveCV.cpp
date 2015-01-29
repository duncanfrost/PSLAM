//Copyright 2008 Isis Innovation Limited 
#include "VideoSourceLiveCV.h"

using namespace std;


VideoSourceLiveCV::VideoSourceLiveCV():amoVideoSource()
{
	cv::VideoCapture cap2;
	cap2.open(0);
	if(!cap2.isOpened())  // check if we succeeded
		std::cerr<<"acquisition problem"<<std::endl;
	
	cv::Mat frame;
	cap2 >> frame;
	
	//mirSize[0]=frame.size().width;
	//mirSize[1]=frame.size().height;	
	mirSize[0]=640;
	mirSize[1]=480;	
	myCamera.Init(mirSize[0],mirSize[1],myCameraType);
	record=false;
	
};
VideoSourceLiveCV::VideoSourceLiveCV(CameraType _camType):amoVideoSource(_camType)
{
	cv::VideoCapture cap2;
	cap2.open(0);
	if(!cap2.isOpened())  // check if we succeeded
		std::cerr<<"acquisition problem"<<std::endl;
	
	cv::Mat frame;
	cap2 >> frame;
	
	if(lvl_pyr_used!=0)
	  	for(int i=0;i<lvl_pyr_used;i++)
			cv::pyrDown( frame,frame, cv::Size( frame.cols/2, frame.rows/2 ) );	

	
	mirSize[0]=frame.size().width;
	mirSize[1]=frame.size().height;	
	myCamera.Init(mirSize[0],mirSize[1],myCameraType);
	record=false;
}

void VideoSourceLiveCV::initCam()
{
  //logically it should be able to do that step at constructor, 
  //however was bugging ... to solve
  cv::VideoCapture *cap;cap=new cv::VideoCapture;
  cap->open(0);
  if(!cap->isOpened())  // check if we succeeded
        std::cerr<<"acquisition problem"<<std::endl;
  cvCap=(void*)cap;
  
  frameId=0;
};

void VideoSourceLiveCV::grabNewFrame()
{
	pthread_mutex_lock( &mutex_image );	
	cv::VideoCapture *cap=(cv::VideoCapture*)cvCap;
	
	cv::Mat frame;
	*cap >> frame;
	
	if(lvl_pyr_used!=0)
	  	for(int i=0;i<lvl_pyr_used;i++)
			cv::pyrDown( frame,frame, cv::Size( frame.cols/2, frame.rows/2 ) );	

	if (frame.cols == 0) {
	    cout << "Error getting live image " << endl;
	}
        cv::resize(frame,imRGB,cv::Size(mirSize[0],mirSize[1]),CV_INTER_LINEAR);
	cvtColor(imRGB,imBW,CV_RGB2GRAY);

	if(record)
		recordFrame();
	pthread_mutex_unlock( &mutex_image );

	frameId++;
}


