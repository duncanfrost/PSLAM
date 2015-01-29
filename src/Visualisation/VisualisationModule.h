#ifndef VISUMODULE_H
#define VISUMODULE_H

#include "../AmoDefines.h"
#include HGLUT

#include "MotherWindow.h"
#include "EmptyWindow.h"
#include "MapWindow.h"
#include "SimuWindow.h"
#include "PlottingWindow.h"
#include "ImageWindow.h"
#include "../ImageSource/VideoSource.h"

#include <list>




class VisualisationModule
{
 public:
	VisualisationModule(void (*_idle_fct)());
	
	//create a window displaying an image, can set processKey function and draw additional stuff
	void addWindowImage(std::string _title,amoVideoSource *_vid_source);
	void addWindowImage(std::string _title,cv::Mat *_img);
	void addWindowImage(std::string _title,cv::Mat *_img,MotherWindow *&_res_window);
	//void addWindowImage(std::string _title,cv::Mat *_img,void (*_processKeys_fct)(unsigned char key, int x, int y),void (*_processMouse_fct)(int button, int state, int x, int y) ,void (*_processMouseActiveMotion_fct)(int x, int y) ,void (*_addDrawFunction)());
	//create an empty window with draw function to define
	void addWindowEmpty(std::string _title,int _width,int _height);
	void addWindowEmpty(std::string _title,int _width,int _height,EmptyWindow *&_res_window);
	//create a window displaying current frame and reprojection of estimated map
	void addWindowMap(std::string _title,int _width,int _height,Camera *_myCamera,MapWindow *&_res_window);
	void addWindowSimu(std::string _title,std::string _objBaseName,std::string _objFile,int _width,int _height,Camera *_myCamera,SimuWindow *&_res_window);
	//create a window to display plot
	void addWindowPlot(std::string _title,int _width,int _height,PlottingWindow *&_res_window,int nb_plots,int _timeLine,bool _incTimeWhenDisplay=true);
	
	//set event function to last added window
	void setOnDraw(void (*_DrawFunction)());
	void setOnKeyPress(void (*_processKeys_fct)(unsigned char key, int x, int y));
	void setOnMouseClick( void (*processMouse_fct)(int button, int state, int x, int y));
	void setOnMouseMotion( void (*processMouseActiveMotion_fct)(int x, int y));
	
	//init GL , start rendering and processing idle function
	void prepareLoop(int argc, char** argv);
	void startLoop(int argc, char** argv);
	void drawWindows();
	
	~VisualisationModule();

 protected: 
	//std::list<amoWindow> listWindows;
	MotherWindow *listWindows[NB_WINDOWS_MAX];
	int nb_windows;
	void (*idle_fct)();
	//set the position of the windows automatically
	bool autoPositioning;
};

#endif
