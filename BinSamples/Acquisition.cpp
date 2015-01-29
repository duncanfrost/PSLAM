//sample to show how to use live images

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 

#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);

//Visualization
VisualisationModule *VisuEngine;
//image acquisition
VideoSourceLiveCV *myVideoSource;
Camera myCamera;//camera object, calibration ...
	

int main(int argc, char** argv)
{
	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->initCam();

	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
	VisuEngine->setOnKeyPress(&processNormalKeys);


	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


void Idle(void) 
{
	myVideoSource->grabNewFrame();
	VisuEngine->drawWindows();
}

void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			break;
	}


	glutPostRedisplay();
}
