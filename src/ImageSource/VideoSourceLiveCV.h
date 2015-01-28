//Copyright 2008 Isis Innovation Limited 
#ifndef VIDEOSOURCELIVECV_H
#define VIDEOSOURCELIVECV_H

#include <objectr3d/VideoSource.h>

class VideoSourceLiveCV : public amoVideoSource
{
 public:
 	//constructor without calibration
	VideoSourceLiveCV();
	//constructor with calibration
	VideoSourceLiveCV(CameraType _camType);  
  
	void initCam();
	void grabNewFrame();
 private:
	void *mptr;
	void *cvCap;
};

#endif