#include "VisualisationModule.h"



VisualisationModule::VisualisationModule(void (*_idle_fct)())
{
	idle_fct=_idle_fct;
	nb_windows=0;
	autoPositioning=true;
}
void VisualisationModule::setOnDraw(void (*_DrawFunction)()){VisuWindowNs::mImageWinStatic[nb_windows-1].addDrawFunctionWindow=_DrawFunction;};	
void VisualisationModule::setOnKeyPress(void (*_processKeys_fct)(unsigned char key, int x, int y)){VisuWindowNs::mImageWinStatic[nb_windows-1].processKeys_fct=_processKeys_fct;};
void VisualisationModule::setOnMouseClick( void (*processMouse_fct)(int button, int state, int x, int y)){VisuWindowNs::mImageWinStatic[nb_windows-1].processMouse_fct=processMouse_fct;};
void VisualisationModule::setOnMouseMotion( void (*processMouseActiveMotion_fct)(int x, int y)){VisuWindowNs::mImageWinStatic[nb_windows-1].processMouseActiveMotion_fct=processMouseActiveMotion_fct;};

void VisualisationModule::addWindowImage(std::string _title,amoVideoSource *_vid_source )
{
	ImageWindow *newWindow=new ImageWindow(_title,_vid_source);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;
}
/*void VisualisationModule::addWindowImage(std::string _title,cv::Mat *_img,MotherWindow *&_res_window,void (*_processKeys_fct)(unsigned char key, int x, int y),void (*_addDrawFunction)() )
{
	//ImageWindow *newWindow=new ImageWindow(_title,_img,_processKeys_fct,_addDrawFunction);
	ImageWindow *newWindow=new ImageWindow(_title,_img,_processKeys_fct,_addDrawFunction);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	_res_window=(MotherWindow*)newWindow;
	nb_windows++;	
}*/
void VisualisationModule::addWindowImage(std::string _title,cv::Mat *_img,MotherWindow *&_res_window)
{
	//ImageWindow *newWindow=new ImageWindow(_title,_img,_processKeys_fct,_addDrawFunction);
	ImageWindow *newWindow=new ImageWindow(_title,_img);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	_res_window=(MotherWindow*)newWindow;
	nb_windows++;	
}
void VisualisationModule::addWindowImage(std::string _title,cv::Mat *_img)
{
	//ImageWindow *newWindow=new ImageWindow(_title,_img,_processKeys_fct,_addDrawFunction);
	ImageWindow *newWindow=new ImageWindow(_title,_img);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;	
}
/*void VisualisationModule::addWindowImage(std::string _title,cv::Mat *_img,void (*_processKeys_fct)(unsigned char key, int x, int y) ,void (*_processMouse_fct)(int button, int state, int x, int y) ,void (*_processMouseActiveMotion_fct)(int x, int y) ,void (*_addDrawFunction)())
{
	ImageWindow *newWindow=new ImageWindow(_title,_img,_processKeys_fct,_processMouse_fct,_processMouseActiveMotion_fct,_addDrawFunction);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;	
}*/
void VisualisationModule::addWindowEmpty(std::string _title,int _width,int _height)
{
	EmptyWindow *newWindow=new EmptyWindow(_title,_width,_height);
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;		
}
void VisualisationModule::addWindowEmpty(std::string _title,int _width,int _height,EmptyWindow *&_res_window)
{
	EmptyWindow *newWindow=new EmptyWindow(_title,_width,_height);
	_res_window=newWindow;
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;		
}


void VisualisationModule::addWindowMap(std::string _title,int _width,int _height,Camera *_myCamera,MapWindow *&_res_window)
{
	MapWindow *newWindow=new MapWindow(_title,_width,_height,_myCamera);
	_res_window=newWindow;
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;		
	
}
void VisualisationModule::addWindowSimu(std::string _title,std::string _objBaseName,std::string _objFile,int _width,int _height,Camera *_myCamera,SimuWindow *&_res_window)
{
	SimuWindow *newWindow=new SimuWindow(_title,_objBaseName,_objFile,_width,_height,_myCamera);
	_res_window=newWindow;
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;		
	
}
void VisualisationModule::addWindowPlot(std::string _title,int _width,int _height,PlottingWindow *&_res_window,int nb_plots,int _timeLine,bool _incTimeWhenDisplay)
{
	PlottingWindow *newWindow=new PlottingWindow(_title,_width,_height,nb_plots,_timeLine,_incTimeWhenDisplay);
	_res_window=newWindow;
	listWindows[nb_windows]=(MotherWindow*)newWindow;
	nb_windows++;			
}

void VisualisationModule::prepareLoop(int argc, char** argv)
{
	glutInit(&argc, argv);	
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	
	glutIdleFunc(*idle_fct);

	
	//InitPostIdle for each window
	//for(int i=0;i<nb_windows;i++)
	for(int i=nb_windows-1;i>=0;i--)
	{
		listWindows[i]->CreateWindow();
		VisuWindowNs::mImageWinStatic[i].windowID=listWindows[i]->get_window();
		listWindows[i]->setEvents();
		//if(VisuWindowNs::mImageWinStatic[nb_windows-1].addDrawFunctionWindow)
			
	}
	
	//position windows in a mosaic
	if(autoPositioning)
	{
		int x_current=0,y_current=0;
		int screenWidth = glutGet(GLUT_SCREEN_WIDTH); 
		for(int i=0;i<nb_windows;i++)
		{
			glutSetWindow(listWindows[i]->get_window());
			glutPositionWindow(x_current,y_current);
			x_current+=listWindows[i]->get_width();
			if(x_current>screenWidth-listWindows[i]->get_width())
			{
				x_current=0;
				y_current+=listWindows[i]->get_height();
			}
		}

	}
}
void VisualisationModule::startLoop(int argc, char** argv)
{
	glutMainLoop();
	
}
void VisualisationModule::drawWindows()
{
	for(int i=0;i<nb_windows;i++)
	{
		glutSetWindow(listWindows[i]->get_window());
		listWindows[i]->prepare_draw();
		listWindows[i]->drawObjects();
		
		if(VisuWindowNs::mImageWinStatic[i].addDrawFunctionWindow)
			(*VisuWindowNs::mImageWinStatic[i].addDrawFunctionWindow)();
		
		glutSwapBuffers();	
		glutPostRedisplay();
	}
}
VisualisationModule::~VisualisationModule()
{
	for(int i=0;i<nb_windows;i++)
	{
		glutDestroyWindow(listWindows[i]->get_window());
		delete listWindows[i];
	}

}
