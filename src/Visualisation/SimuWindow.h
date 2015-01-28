#pragma once

#include <vector>


#include <objectr3d/Camera.h>
#include <objectr3d/MotherWindow.h>
#include <objectr3d/HomogeneousMatrix.h>
#include <objectr3d/obMap.h>
#include <objectr3d/objloader.h>

class SimuWindow:public MotherWindow
{
 public:
	SimuWindow(std::string _title,std::string _objBaseName,std::string _objFile,int _width,int _height,Camera *_myCamera);
	void CreateWindow();
	void prepare_draw();
	void setEvents();
	
	//set the position of the camera from which map is viewed
	void setCameraPose(HomogeneousMatrix _cHc2);
	void moveCamera(HomogeneousMatrix _cHc2);
	
	//set velocity of camera for keyboard control
	void set_velocity_translation(float _f);
	
	void setLight();
	
	//get current data
	HomogeneousMatrix getCurrentPose(){return displayedPose;};
	cv::Mat &getCurrentImage(){return current_img;};
	cv::Mat &getCurrentImageBW(){return current_img_BW;};
	cv::Mat &getCurrentDepth(){return current_img_Z;};
	
	
 protected: 
	 //pointer to camera calibration
	 Camera *myCamera;
	 bool show_texture;
	 
	 std::string objBaseName;
	 std::string objFile;
	 int objId;
	 objloader obj;
	 
	 //current output data
	 HomogeneousMatrix displayedPose;
	 cv::Mat current_img;
	 cv::Mat current_img_BW;
	 cv::Mat current_img_Z;
};
