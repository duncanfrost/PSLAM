//Copyright 2008 Isis Innovation Limited 
#include "VideoSource.h"


amoVideoSource::amoVideoSource(CameraType _camType)
{
	//force resolution to be 640x480
	mirSize[0]=640;
	mirSize[1]=480;

	myCameraType=_camType;
	//myCamera.Init(mirSize[0],mirSize[1],myCameraType);
		
	frameId=0;	
	
	imBW.create(mirSize[1], mirSize[0], CV_8UC1);
	imRGB.create(mirSize[1], mirSize[0], CV_8UC3);
	
	pthread_mutex_init(&mutex_image,NULL);
	lvl_pyr_used=0;
	record=false;
}

void amoVideoSource::setLvlAcquisition(short _l)
{
	lvl_pyr_used=_l;
	if(lvl_pyr_used!=0)
	  	for(int i=0;i<lvl_pyr_used;i++)
			cv::pyrDown( imRGB,imRGB, cv::Size( imRGB.cols/2, imRGB.rows/2 ) );	
	cv::cvtColor(imRGB,imBW,CV_RGB2GRAY);
	mirSize[0]=imBW.size().width;
	mirSize[1]=imBW.size().height;
	myCamera.Init(mirSize[0],mirSize[1],myCameraType);
	
}

void amoVideoSource::GetBWFrame(cv::Mat &_imBW)  
{
	//_imBW.create(mirSize[1], mirSize[0], CV_8UC1, cv::Scalar(0));
	pthread_mutex_lock( &mutex_image );	
	_imBW=imBW.clone();
	pthread_mutex_unlock( &mutex_image );
}

void amoVideoSource::GetRGBFrame(cv::Mat &_imRGB) 
{
	//_imRGBA.create(mirSize[1], mirSize[0], CV_8UC3, cv::Scalar(0));
	pthread_mutex_lock( &mutex_image );	
	_imRGB=imRGB.clone();
	pthread_mutex_unlock( &mutex_image );
}
  
void amoVideoSource::recordFrame() 
{
	char fileName[200];
	sprintf(fileName, printfPathRecord.c_str(), frameId);
	cv::imwrite(fileName,imRGB);
}
  
amoVideoSource::~amoVideoSource()
{

}

