//Copyright 2008 Isis Innovation Limited 
#ifndef IMAGEWINDOW_H
#define IMAGEWINDOW_H


#include "MotherWindow.h"
#include "glTexture.h"
#include "../ImageSource/VideoSource.h"

//TextureSet imgTexture;
	
class ImageWindow:public MotherWindow
{
 public:
	ImageWindow(std::string _title,amoVideoSource *_vid_source);
	ImageWindow(std::string _title,cv::Mat *_img);
	void CreateWindow();
	void prepare_draw();
	
	void InitPreIdle(){};
	void InitPostIdle(){};

 protected: 
	amoVideoSource *vid_source;
	cv::Mat *pt_img;
};

#endif
