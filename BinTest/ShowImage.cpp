//sample for use of GUI

#define AMOVERBOSE 1

#include <fstream>
#include <objectr3d/Camera.h>

#include <objectr3d/VisualisationModule.h>

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);
void processNormalKeys2(unsigned char key, int x, int y);
void addDrawFunction(void) ;

//Visualization
VisualisationModule *VisuEngine;

void onClick(void) ;
	
cv::Mat Img1;	
cv::Mat Img2;	
	
MotherWindow *window1;
MotherWindow *window2;

int main(int argc, char** argv)
{
	Img1=cv::imread("/home/amaury/Pictures/jump.png");
	Img2=cv::imread("/home/amaury/Pictures/jump2.png");

	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("image1",&Img1,window1);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeys);
	VisuEngine->addWindowImage("image2",&Img2,window2);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeys2);

	VisuEngine->prepareLoop(argc, argv);
	
	window1->addButton("Button",&onClick);
	
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}

void onClick()
{
	std::cout<<"click"<<std::endl;
}

void Idle(void) 
{
	VisuEngine->drawWindows();
}

void processNormalKeys(unsigned char key, int x, int y)
{
	std::cout<<"processNormalKeys"<<std::endl;;
	switch(key) {
		case 27://esc
			exit(0);
			break;
	}


	glutPostRedisplay();
}
void processNormalKeys2(unsigned char key, int x, int y)
{
	std::cout<<"processNormalKeys2"<<std::endl;;
	switch(key) {
		case 27://esc
			std::cout<<"esc"<<std::endl;;
			break;
	}


	glutPostRedisplay();
}

void addDrawFunction(void) 
{	
 
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glEnable(GL_POINT_SMOOTH);
	
	set2DGLProjection();
	glColor3f(0,1,0);
	

				glBegin(GL_LINES);	
				glVertex2f(0,0);
				glVertex2f(100,100);
				glEnd();				

	glColor3f(1,1,1);
	unset2DGLProjection();

}