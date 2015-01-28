//estimate rotation between two consecutive frames

#define AMOVERBOSE 1

#include <fstream>
#include <objectr3d/Camera.h>

#include <objectr3d/VisualisationModule.h>
#include <objectr3d/VideoSourceLiveCV.h>
#include <objectr3d/CvWrapper.h>
#include <objectr3d/ImageProcess.h>
#include <objectr3d/GlobalTransfoEstim.h>

#include <Eigen/Core>
using namespace Eigen;

int main(int argc, char** argv)
{
	InitProcAndGPU();
	
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






