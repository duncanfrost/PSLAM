/*============================================================================
  RotationEstim.cpp: Sample program that estimates the rotation between two
  images. Not currently working.
============================================================================*/

#define AMOVERBOSE 1

#include <fstream>


#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"
#include "../src/CvWrapper.h"
#include "../src/TrackEngines/ImageProcess.h"
#include "../src/TrackEngines/GlobalTransfoEstim.h"

#include <Eigen/Core>
using namespace Eigen;

int main(int argc, char** argv)
{
    std::cout << "Got here" << std::endl;
	InitProcAndGPU();
    std::cout << "Got here" << std::endl;
	
	cv::Mat Img0=cv::imread("/home/amaury/Pictures/jump.png",CV_LOAD_IMAGE_GRAYSCALE);
	if (Img0.cols == 0) {
	    std::cout << "Error reading image " << std::endl;
	}
	
	cv::Mat Img1;GetGauss(Img0,Img1);
	cv::Mat Img;GetGauss(Img0,Img);//320x240

	Camera myCamera(Img.size().width,Img.size().height,CamPlaystationEye);

	//rotate image
	Vector4f p;p[0]=-15;p[1]=0;p[2]=0.2;p[3]=0.1;
	//Vector4f p;p[0]=0;p[1]=0;p[2]=0.;p[3]=0.05;
	std::cout<<"GT tranfo = "<<p.transpose()<<std::endl;
	cv::Mat ImgRot=applySE2TzImage(Img,p,&myCamera);
	cv::imwrite("rot.jpg",ImgRot);
	
	amoTimer timer;timer.start();
	Vector4f pEst;
	float residue=estimateSE2TzPyr(Img,ImgRot,&myCamera,pEst,3,2,6);
	std::cout<<"Est tranfo = "<<pEst.transpose()<<std::endl;
	
	cv::Mat Imgtest=applySE2TzImage(Img,pEst,&myCamera);
	cv::imwrite("rotEst.jpg",Imgtest);
	

	HomogeneousMatrix Hest=SE2TztoSO3(pEst,1,&myCamera);
	std::cout<<"estimed transfo = "<<Hest.get_p().transpose()<<std::endl;
	timer.stop();


	return 0;
}






