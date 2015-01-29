//there are three kinds of classes originally in libViso
//this sample shows the extracted features and there class

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);
void addDrawFunction(void) ;
void useNewMatches(std::vector<p_match> &_matchesCurrent,std::vector<p_match> &_matchesCurrentBucket);
void EstimateCurrentPose(std::vector<p_match> &_matchesCurrent,HomogeneousMatrix22 &relPoseBA);

//Visualization
VisualisationModule *VisuEngine;
MapWindow *mMapWindow;
//image acquisition
VideoSourceLiveCV *myVideoSource;
Camera myCamera;//camera object, calibration ...

//reference
HomogeneousMatrix22 PoseRef;
cv::Mat imgRef;	

cv::Mat imgDisplay;	

//tracking info
Matcher matcher;
void initMatcher();



std::vector<p_feat> feature_p[4];

int main(int argc, char** argv)
{
	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->initCam();
	myCamera=*myVideoSource->getPointerCamera();

	//get reference image
	myVideoSource->grabNewFrame();
	myVideoSource->GetFramePointerBW()->copyTo(imgRef);
	initMatcher();
	
	//init inv depth features
	for(int i=0;i<4;i++)
		feature_p[i]=matcher.getFeaturespClass(i);
	
	//alloc AllMatches
	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",&imgRef);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeys);


	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}

bool poseProcess=false;
void Idle(void) 
{
	VisuEngine->drawWindows();
}


void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			break;
		case ' '://esc
			poseProcess=!poseProcess;
			break;
	}


	glutPostRedisplay();
}

void initMatcher()
{
	
	uint8_t* img_data  = imgRef.data;
	
	// image dimensions
	int32_t dims[3];
	dims[0] = imgRef.cols; // width
	dims[1] = imgRef.rows; // height
	dims[2] = dims[0]; // image width, BW so it equals to width
			
	matcher.pushBack(img_data,dims,true);	
	
	//coutRed<<"Init Matcher"<<endlRed;
	//std::vector<p_feat> feature_current=matcher.getFeatures(img_data,dims);
}

void addDrawFunction(void) 
{	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	set2DGLProjection();
	glPointSize(3.0);
	glColor3f(0,1,0);

	for(int c=0;c<4;c++)	
	{
		if(c==0)glColor3f(1,0.,0);
		if(c==1)glColor3f(0,1.,0);
		if(c==2)glColor3f(0,0.,1);
		if(c==3)glColor3f(0,1.,1);
		for(int i=0;i<feature_p[c].size();i++)
		{
			glLineWidth(1.0);
			glBegin(GL_POINTS);	
			glVertex2f(feature_p[c][i].u1c,feature_p[c][i].v1c);
			glEnd();
		}
	}

	
	glColor3f(1,1,1);
	unset2DGLProjection();

}
