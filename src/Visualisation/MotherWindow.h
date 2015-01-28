#ifndef AMOWINDOW_H
#define AMOWINDOW_H

#include <string>
#include <objectr3d/AmoDefines.h>
#include HGLUT
#include <objectr3d/glTools.h>
#include <objectr3d/VisuButton.h>

#define NB_WINDOWS_MAX 10

//buttons in windows

//list of event functions that can be added to normal windows events
namespace VisuWindowNs
{
	struct VisuWindowStatic
	{
		int windowID;
		void (*addDrawFunctionWindow)();
		void (*processKeys_fct)(unsigned char key, int x, int y);
		void (*processMouse_fct)(int button, int state, int x, int y);
		void (*processMouseActiveMotion_fct)(int x, int y);
		VisuWindowStatic()
		{
		  windowID=-1;
		  addDrawFunctionWindow=NULL;
		  processKeys_fct=NULL;
		  processMouse_fct=NULL;
		  processMouseActiveMotion_fct=NULL;
		}
	};

	extern std::vector<VisuButton> mButtons[NB_WINDOWS_MAX];
	extern VisuWindowStatic mImageWinStatic[NB_WINDOWS_MAX];
	int windowIDtoStatic(int _windowID);
}

class MotherWindow
{
 public:
	MotherWindow();
	~MotherWindow(){};
	
	virtual void CreateWindow()=0;
	virtual void setEvents();
	
	//function that is called before each drawing, usefull since draw is static while prepare_draw is not
	virtual void prepare_draw()=0;
	void drawObjects();
	
	int get_window(){return window;};
	int get_width(){return width;};
	int get_height(){return height;};
	
	void setTitle(const char *_title);
	
	//add button?
	void addButton(std::string _caption,void (*_onClick)());
	
	//mouse control in case of button
	

 protected: 
	std::string title;
	int window;
	int pos_x,pos_y;
	int width,height;

};

#endif