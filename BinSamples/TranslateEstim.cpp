//estimate translation between two consecutive frames

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
	//get one camera calibration, necessary for back projection of image plane
	Camera myCamera(640,480,CamPlaystationEye);
	
	cv::Mat Img=cv::imread("/home/adame/im2.jpg",CV_LOAD_IMAGE_GRAYSCALE);
	cv::Mat Img2;GetGauss(Img,Img2);
	
	//translate image
	Vector2f t;t[0]=20;t[1]=5;
	cv::Mat Imgt2=translateImage(Img2,t);
	cv::imwrite("rot.jpg",Imgt2);
	
	Vector2f tEst=estimateTranslation(Img2,Imgt2,4);
	std::cout<<"trans = "<<tEst.transpose()<<std::endl;

	cv::Mat Imgtest=translateImage(Img2,tEst);
	cv::imwrite("rotEst.jpg",Imgtest);

	return 0;
}






