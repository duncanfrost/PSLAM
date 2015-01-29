//acquire live images, save it on drive when space bar is pressed
//could have use myVideoSource->setRecord("/media/Data/Datas/NewAcquisition/img_%d.png");
//but here we put images in cash to have better frame rate


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
	
#define MAX_TRAJ 5000
#define DIRNAME "/media/Data/Datas/NewAcquisition/"
cv::Mat buf_img[MAX_TRAJ];

int main(int argc, char** argv)
{
	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->initCam();
	

	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
	VisuEngine->setOnKeyPress(&processNormalKeys);

	std::cout<<"Press Space bar to stop and save sequence to disk."<<std::endl;
	
	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}

int cpt_img=0;

void Idle(void) 
{
	myVideoSource->grabNewFrame();
	VisuEngine->drawWindows();
	
	if(cpt_img<MAX_TRAJ)
	{
		myVideoSource->GetRGBFrame(buf_img[cpt_img]);
		cpt_img++;
	}
}

void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			break;
		case ' '://esc
			std::cout<<"Storing to disk"<<std::endl;
			for(int i=0;i<cpt_img;i++)
			{
				if(i%(cpt_img/10)==0)std::cout<<100.*i/cpt_img<<" %"<<std::endl;
				char filename_img[200];
				sprintf (filename_img, "%simg_%d.png",DIRNAME, i);
				cv::imwrite(filename_img,buf_img[i]);				
			}

			exit(0);
			break;
			
	}


	glutPostRedisplay();
}
