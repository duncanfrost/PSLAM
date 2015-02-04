//Copyright 2008 Isis Innovation Limited 
#ifndef VIDEOSOURCE_H
#define VIDEOSOURCE_H


#include "../AmoDefines.h"
#include "../Primitives/Camera.h"

#include <pthread.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>



//define projection inverse: get 3D point world coordinate from its pixel coord, depth and camera modelview (world coord to cam coord matrix)
class amoVideoSource;


class amoVideoSource
{
 public:
	//constructor without calibration
	//amoVideoSource();
	//constructor with calibration
	amoVideoSource(CameraType _camType=CamUndefined);
	void initCam();
	
	//acquire a new frame
	virtual void grabNewFrame()=0;
	
	//get the frame previously acquired (with copy)
	void GetBWFrame(cv::Mat &_imBW);  
	void GetRGBFrame(cv::Mat &_imRGB); 
	
	//get pointer to the frame previously acquired (no copy)
    cv::Mat* GetFramePointer(){return &imRGB;}
    cv::Mat* GetFramePointerBW(){return &imBW;}
	
	//get image size
    Vector2i getSize(){return mirSize;}
    int getWidth(){return mirSize[0];}
    int getHeight(){return mirSize[1];}
	//get frame number
    int getFrameId(){return frameId;}

    Camera* getPointerCamera(){return &myCamera;}
	~amoVideoSource();

	//if want to use lvl pyramid instead of raw image
	void setLvlAcquisition(short _l);
	
    void setRecord(std::string _printfPathRecord){printfPathRecord=_printfPathRecord;record=true;}
	void recordFrame();

 protected: 
	 //acquisition props
	int frameId;
	Vector2i mirSize;
	
	//acquired image from GetFrame, return by GetFramePointer
	cv::Mat imBW ;
	cv::Mat imRGB ;

	//GetFrame and grabNewFrame might be in different threads so access to images is protected by mutex
	pthread_mutex_t mutex_image;		

	//camera props
	Camera myCamera;
	CameraType myCameraType;
	
	short lvl_pyr_used;
	
	bool record;
	std::string printfPathRecord;
};

#endif
