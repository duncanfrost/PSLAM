#pragma once

#include <vector>


#include "../Primitives/Camera.h"
#include "MotherWindow.h"
#include "../Primitives/HomogeneousMatrix.h"
#include "../Primitives/obMap.h"

#define SHOW_BEST_KF_PAIR
//#define SHOW_BAD_POINTS

//TextureSet imgTexture;
struct PointInvDepth
{
    int srcKf;
    Vector2f meterCoord;
    float invDepth;
    float invDepthCovar;  
    PointInvDepth(){invDepthCovar=0;}
};
	

class MapWindow:public MotherWindow
{
 public:
	MapWindow(std::string _title,int _width,int _height,Camera *_myCamera);
	void CreateWindow();
	void prepare_draw();
	void setEvents();
	
	void addCamera(HomogeneousMatrix *_camPose,Vector3f _col,float _lineSize=1.);
	void addPointCloud(std::vector<Vector3f> *_map);
    void addPointCloud(std::vector<PointInvDepth> *_map){pointInvClouds.push_back(_map);}
	void addMap(obMap *_map);//if several maps are added then use first one as GT
	
	//set the position of the camera from which map is viewed
	void setCameraPose(HomogeneousMatrix _cHc2);
	void moveCamera(HomogeneousMatrix _cHc2);
	HomogeneousMatrix getCameraPose();
	
	//set velocity of camera for keyboard control
	void set_velocity_translation(float _f);
	void set_camera_drawn_size(float _f);
	
    void showClosestKF(int i){closestKF=i;}
    void showActiveKF(std::vector<int> &_activeKF){activeKF=_activeKF;}
	
    void showTextures(){b_showTextures=!b_showTextures;}
    void showFeatureConnections(){b_showFeatureConnections=!b_showFeatureConnections;}
    void showLocalFeature(){b_showLocalFeature=!b_showLocalFeature;}


	
 protected: 
	 //pointer to camera calibration
	 Camera *myCamera;
	 
	//list of cameras poses to display
	std::vector<HomogeneousMatrix*> cameraPoses;
	std::vector<Vector3f> colCamera;
	std::vector<float> lineSizeCamera;
	
	//list of point clouds to display
	std::vector<std::vector<Vector3f>* > pointClouds;
	std::vector<std::vector<PointInvDepth>* > pointInvClouds;
	
	//list of complete map (regroups point cloud and camera keyFrames)
	std::vector<Vector3f> colMapsNew;
	std::vector<obMap*> MapsNew;
	
	std::vector<Vector3f> lotsOfRandColors;

	//shows tracker info
	int closestKF;
	std::vector<int> activeKF;
	
	bool b_showTextures;
	bool b_showFeatureConnections;
	bool b_showLocalFeature;

};
