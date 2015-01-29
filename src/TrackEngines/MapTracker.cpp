
#include "MapTracker.h"
//#include "../VisuOdo/motionEstim.h"



MapTracker::MapTracker(Camera *_cam,obMap *_map)
{

    std::cout << "Got here" << std::endl;
    std::cout << "Size: " << id_closestKF.size() << std::endl;
    std::cout << "Address: " << &id_closestKF << std::endl;
    std::cout << "Address of maptracker: " << this << std::endl;


    coutGreen << "=====================" << endlGreen;
    coutGreen << sizeof(relPose) << endlGreen;
    coutGreen << sizeof(HomogeneousMatrix22) << endlGreen;
    coutGreen << sizeof(Matrix4f) << endlGreen;
    coutGreen << sizeof(float) << endlGreen;
    coutGreen << "=====================" << endlGreen;

    MotionPriorHomography(0,0) = 0.0f;
    //relPose = new HomogeneousMatrix22;
    Init(_cam,_map);
    Test();
}

void MapTracker::Init(Camera *_cam,obMap *_map)
{
	myCamera=_cam;
	myMap=_map;	
	idrelKF=-1;
    MotionPriorHomography(0,0) = 1;
    MotionPriorHomography(1,2) = 1;
    MotionPriorHomography(2,1) = 1;
    std::cout << "Address: " << &id_closestKF << std::endl;

}


