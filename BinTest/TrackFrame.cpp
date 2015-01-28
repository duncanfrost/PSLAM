//define first acquired image as reference and do feature matching
//using libViso with current image and estimate homography between 
//two images, warp current image on ref

#define AMOVERBOSE 1

#include <fstream>
#include <objectr3d/Camera.h>

#include <objectr3d/VisualisationModule.h>
#include <objectr3d/VideoSourceLiveCV.h>
#include <objectr3d/matcherWrapper.h>

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);
void addDrawFunction(void) ;

//Visualization
VisualisationModule *VisuEngine;
//image acquisition
VideoSourceLiveCV *myVideoSource;
Camera myCamera;//camera object, calibration ...

//reference
#define nb_lvl_img 3
cv::Mat imgRef[nb_lvl_img];	
cv::Mat imgCurr[nb_lvl_img];	
cv::Mat imgCurr0_warp;	

//tracking info viso
matcherWrapper matcher;
Matrix3f Homography;//current homography with ref
std::vector<p_match> matchesCurrent;

int main(int argc, char** argv)
{
	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->initCam();
	myVideoSource->setLvlAcquisition(1);

	//get reference image
	myVideoSource->grabNewFrame();
	myVideoSource->GetFramePointerBW()->copyTo(imgRef[0]);
	for(int i=1;i<nb_lvl_img;i++)
		cv::pyrDown( imgRef[i-1], imgRef[i], cv::Size( imgRef[i-1].cols/2, imgRef[i-1].rows/2 ) );	

	imgRef[0].copyTo(imgCurr0_warp);
	
	matcher.InitWithRef(imgRef);
	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeys);
	VisuEngine->addWindowImage("Input image warped",&imgCurr0_warp);
	VisuEngine->setOnKeyPress(&processNormalKeys);


	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}

bool pauseProc=false;
void Idle(void) 
{
	if(!pauseProc)
	{
		std::cout<<"New image"<<std::endl;
		myVideoSource->grabNewFrame();
		myVideoSource->GetFramePointerBW()->copyTo(imgCurr[0]);
		for(int i=1;i<nb_lvl_img;i++)
			cv::pyrDown( imgCurr[i-1], imgCurr[i], cv::Size( imgCurr[i-1].cols/2, imgCurr[i-1].rows/2 ) );	
		

		amoTimer timerIdle;
		timerIdle.start();
		
		matcher.match(imgCurr);
		matchesCurrent=matcher.getCurrentMatches();
		Homography=matcher.getHomography();
		std::cout<<"matchesCurrent.size() = "<<matchesCurrent.size()<<std::endl;
		timerIdle.stop("Idle");
		
		imgCurr0_warp= warpImageInv(imgCurr[0],Homography);
	
		VisuEngine->drawWindows();
	}
}

void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			break;
		case ' ':
			pauseProc=!pauseProc;
			break;
	}


	glutPostRedisplay();
}


void addDrawFunction(void) 
{	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	set2DGLProjection();
	glPointSize(3.0);
	glLineWidth(1);
	glColor3f(0,1,0);


	for(int i=0;i<matchesCurrent.size();i++)
	{
		glLineWidth(1.0);
		glColor3f(0,0.,1);
		glBegin(GL_LINES);	
		glVertex2f(matchesCurrent[i].u1p,matchesCurrent[i].v1p);
		glColor3f(0,1.,1);
		glVertex2f(matchesCurrent[i].u1c,matchesCurrent[i].v1c);
		glEnd();
	}

	
	glColor3f(1,1,1);
	unset2DGLProjection();

}
