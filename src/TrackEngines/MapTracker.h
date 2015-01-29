//Map tracking module.

#pragma once

//#include <opencv/cv.h>
//#include <opencv/highgui.h>
//#include <sys/time.h>

#include <Eigen/Core>
using namespace Eigen;

#include "../AmoDefines.h"
#include "../Primitives/Camera.h"
//#include "ImageProcess.h"
#include "../Primitives/obMap.h"
#include "../Primitives/HomogeneousMatrix.h"
//#include "GlobalTransfoEstim.h"
//#include "../MapEngines/MapOptimiser.h"
//#include "../MapEngines/MapOptimiserEssential.h"

//for fast corner detection
//#include <opencv/cvaux.h> //for opencv 2.4
//#include <opencv2/opencv.hpp> // for later one


class MapTracker
{
public:
    MapTracker()
    {

    };
	~MapTracker(){};
    MapTracker(Camera *_cam,obMap *_map);

	void Init(Camera *_cam,obMap *_map);
    void Test(void)
    {
        coutGreen << "=====================" << endlGreen;
        coutGreen << sizeof(relPose) << endlGreen;
        coutGreen << sizeof(HomogeneousMatrix22) << endlGreen;
        coutGreen << sizeof(Matrix4f) << endlGreen;
        coutGreen << sizeof(float) << endlGreen;
        coutGreen << "=====================" << endlGreen;
    }
	




	int getIdRelKF(){return idrelKF;};
    std::vector<int> getClosestKFs()
    {

        coutGreen << "=====================" << endlGreen;
        coutGreen << sizeof(relPose) << endlGreen;
        coutGreen << sizeof(HomogeneousMatrix22) << endlGreen;
        coutGreen << sizeof(Matrix4f) << endlGreen;
        coutGreen << sizeof(float) << endlGreen;
        coutGreen << "=====================" << endlGreen;



        std::cout << "Address of maptracker: " << this << std::endl;

        std::vector<int> out;
        std::cout << "Size: " << id_closestKF.size() << std::endl;
        std::cout << "Address: " << &id_closestKF << std::endl;

        std::cout << sizeof(relPose) << std::endl;
         std::cout << sizeof(HomogeneousMatrix22) << std::endl;





        for (unsigned int i = 0; i < id_closestKF.size(); i++)
        {
            int t = id_closestKF[i];
            out.push_back(t);
        }

        return out;
    }
private:
  
	//pointers o camera and map
	Camera *myCamera;
	obMap *myMap;
	
	//current relative pose between KFcurrent and KFclosest (in map)
    HomogeneousMatrix22 relPose;
	
	int idrelKF;//id of KF relative pose is defined
	std::vector<int> id_closestKF;//id of current closest KF in map (active window)
	
	Matrix3f MotionPriorHomography;//can use MotionPriorHomography=Id or use constant mvt prior=> MotionPriorHomography=Homography(t-1)*Homography(t-2).inverse()
	//for instance in Kitti, when rotate, rotation is quiet smooth but large=> if we don not take motion prior into account lots of features are not matched
	
};

	
	
