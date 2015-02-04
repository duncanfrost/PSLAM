//Copyright 2008 Isis Innovation Limited 
#ifndef VIDEOSOURCESEQ_H
#define VIDEOSOURCESEQ_H

#include "VideoSource.h"

class VideoSourceSeq : public amoVideoSource
{
 public:
 	//constructor without calibration
	VideoSourceSeq();
	//constructor with calibration
	VideoSourceSeq(char *_printfPath,CameraType _camType,int id0=0);  
  
	void initCam();
	void grabNewFrame();
	
    void setIncrement(int _i){increment=_i;}
 private:
	char *printfPath;
	bool end_sequence;
	
	int increment;
};

#endif
