//Copyright 2008 Isis Innovation Limited 
#include "VideoSourceSeq.h"

using namespace std;


VideoSourceSeq::VideoSourceSeq():amoVideoSource()
{
	frameId=0;
	end_sequence=false;
	printfPath="";
	mirSize=Vector2i(0,0);
	increment=1;
	record=false;

};
//example of path :"/dir/img_%04d.png"
VideoSourceSeq::VideoSourceSeq(char *_printfPath,CameraType _camType,int id0):amoVideoSource(_camType)
{
	frameId=id0;
	printfPath=_printfPath;
	end_sequence=false;
	
	//open first image to get size
	char fileName[200];
	sprintf(fileName, printfPath, frameId);
	imRGB=cv::imread(fileName);
	
	if(lvl_pyr_used!=0)
	  	for(int i=0;i<lvl_pyr_used;i++)
			cv::pyrDown( imRGB,imRGB, cv::Size( imRGB.cols/2, imRGB.rows/2 ) );	

	//imBW=cv::imread(fileName,CV_LOAD_IMAGE_GRAYSCALE);
	cv::cvtColor(imRGB,imBW,CV_RGB2GRAY);
	mirSize[0]=imBW.size().width;
	mirSize[1]=imBW.size().height;
	myCamera.Init(mirSize[0],mirSize[1],myCameraType);
	increment=1;
	record=false;
}

void VideoSourceSeq::initCam()
{

};

#include <iostream>
#include <fstream>
void VideoSourceSeq::grabNewFrame()
{
	if(!end_sequence)
	{
		frameId+=increment;
		char fileName[200];
		sprintf(fileName, printfPath, frameId);

		std::ifstream fout;
		fout.open(fileName);
		if(!fout.is_open())
		{
			std::cerr<<"VideoSourceSeq : End sequence"<<std::endl;
			end_sequence=true;
		}
		fout.close();		
		
		if(!end_sequence)
		{
			imRGB=cv::imread(fileName);
			
			if(lvl_pyr_used!=0)
				for(int i=0;i<lvl_pyr_used;i++)
					cv::pyrDown( imRGB,imRGB, cv::Size( imRGB.cols/2, imRGB.rows/2 ) );	
			cv::cvtColor(imRGB,imBW,CV_RGB2GRAY);
			
			/*imBW.create(mirSize[1], mirSize[0], CV_8UC1);
			for(int y=0;y<mirSize[1];y++)
			for(int x=0;x<mirSize[0];x++)
				imBW.at<uchar>(y, x)=((float)imRGB.at<cv::Vec3b>(y, x)[0]+(float)imRGB.at<cv::Vec3b>(y, x)[1]+(float)imRGB.at<cv::Vec3b>(y, x)[2])/3.;*/
			pthread_mutex_unlock( &mutex_image );
		}
		//if(record)
		//	recordFrame();
	}
}


