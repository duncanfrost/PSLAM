#include "MotherWindow.h"

namespace VisuWindowNs
{
std::vector<VisuButton> mButtons[NB_WINDOWS_MAX];
VisuWindowStatic mImageWinStatic[NB_WINDOWS_MAX];
}

int VisuWindowNs::windowIDtoStatic(int _windowID)
{
	for(int i=0;i<NB_WINDOWS_MAX;i++)
	{
		if(mImageWinStatic[i].windowID==_windowID)
			return i;
	}
	return -1;  
}

using namespace VisuWindowNs;

MotherWindow::MotherWindow()
{

}

void MotherWindow::setTitle(const char *_title)
{
	glutSetWindow(window);
	glutSetWindowTitle(_title);
}

#define default_button_width 90
#define default_button_height 40
#define margin_button 10

void MotherWindow::addButton(std::string _caption,void (*_onClick)())
{
	//put them on top right of window - margin
	int idS=VisuWindowNs::windowIDtoStatic(window);
	
	VisuButton newButton(_caption,_onClick,width-default_button_width-margin_button,default_button_height*mButtons[idS].size()+(1+mButtons[idS].size())*margin_button,default_button_width,default_button_height);
	mButtons[idS].push_back(newButton);
		
}

void MotherWindow::drawObjects()
{
	
	int idS=VisuWindowNs::windowIDtoStatic(window);
	//std::cout<<"Draw object win "<<idS<<" nb buttons = "<<mButtons[idS].size()<<std::endl;
	//std::cout<<"Draw object win "<<idS<<" nb buttons = "<<VisuWindowNs::mImageWinStatic[idS].nb_buttons<<std::endl;
	for(int b=0;b<mButtons[idS].size();b++)
		mButtons[idS][b].draw();
}


void mouse_click(int button, int state, int x, int y)
{
	int window=glutGetWindow();
	int idS=VisuWindowNs::windowIDtoStatic(window);
	
	//do normal thing to do in window => check button action
	for(int b=0;b<mButtons[idS].size();b++)
		mButtons[idS][b].checkForAtivity(button,state,x,y);
	
	//do additional stuff defined in main
	if(VisuWindowNs::mImageWinStatic[idS].processMouse_fct)
		(*VisuWindowNs::mImageWinStatic[idS].processMouse_fct)(button,state,x,y);

}

void mouse_motion(int x, int y)
{
    int window=glutGetWindow();
    int idS=VisuWindowNs::windowIDtoStatic(window);
    
    //do normal thing to do in window => check button action
 	for(int b=0;b<mButtons[idS].size();b++)
		mButtons[idS][b].checkForAtivityMotion(x,y);
   
    //do additional stuff defined in main
    if(VisuWindowNs::mImageWinStatic[idS].processMouseActiveMotion_fct)
	    (*VisuWindowNs::mImageWinStatic[idS].processMouseActiveMotion_fct)(x,y);
}


void MotherWindow::setEvents()
{
	//get id in static array
	int idS=windowIDtoStatic(window);
	if(mImageWinStatic[idS].processKeys_fct)
		glutKeyboardFunc(mImageWinStatic[idS].processKeys_fct);
 	/*if(mImageWinStatic[idS].processMouse_fct)
		glutMouseFunc(mImageWinStatic[idS].processMouse_fct);
	if(mImageWinStatic[idS].processMouseActiveMotion_fct)
		glutMotionFunc(mImageWinStatic[idS].processMouseActiveMotion_fct);*/
	
	glutMouseFunc(mouse_click);
	glutMotionFunc(mouse_motion);
     
}