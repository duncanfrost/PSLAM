//camera object
//load calibration with PTAM style file

//TODO : needs to add radial distortion

#pragma once

#include "../AmoDefines.h"

#include <Eigen/Core>

enum CameraType {CamPlaystationEye,
		  KittiCam,
		  CamFlea2,
		  CamUndefined
};

using namespace Eigen;

class Camera
{
public:
	Camera(){};
	Camera(int _width, int _height, CameraType _camType=CamUndefined);
	void Init(int _width, int _height, CameraType _camType=CamUndefined);
	void InitParams(int _width, int _height);
    ~Camera(){}

	//get calib of camera corresponding to half size image
	Camera getHalfSizeCamera();

	
	// from meter coord in camera frame or z=1 plane to pixels	
	Vector2f Project(const Vector2f& mvLastDistCam);
	Vector2f Project(const Vector3f& mvLastDistCam);
    Vector2f ToPixels(const Vector2f& mvLastDistCam){return Project(mvLastDistCam);}
	Vector2f ToMeters(const Vector2f& mvLastDistCam);
	//project 3D point in Z=1 plane in meters
	Vector2f ProjectZ1(const Vector3f& mvLastDistCam);
	//jacobian of projection (in meter) with respect to 3D point being projected
	MatrixXf ProjectZ1_Jac_X(const Vector3f& mvLastDistCam);
	//jacobian of projection with respect to Delta p applied to camera pose
	MatrixXf ProjectZ1_Jac_Dp(const Vector3f& mvLastDistCam);
	// Inverse operation: pixel to meter
	Vector2f UnProject(const Vector2f& imframe); 
	Vector3f UnProjectZ1(const Vector2f& imframe); // Inverse operation on z=1 plane
	//simply add z=1 to 2d float vector to return its homogeneous representation
    Vector3f ToHomogeneous(const Vector2f& imframe){return Vector3f(imframe[0],imframe[1],1);}
	
	Matrix2f m2PixProjJac();//matrix to convert Jacobian of proj of point in meter to Jacobian in pixel units
	Matrix2f Pix2mProjJac();
    Vector2f getCenterImage(){return mvCenter;}
    Vector2f getFocalPix(){return mvFocal;}
	
	Matrix4f getGlProjectionMat(double near, double far);
    float *getParameters(){return &params[0];}
    float *getParametersPTAM(){return mgvvCameraParams;}
	
	//get approximate size of one pixel in meter
    float pixMeterSize(){return 1./mvFocal[0];}

    int get_width(){return mvImageSize[0];}
    int get_height(){return mvImageSize[1];}
private:
	//parameters from calib file
	float mgvvCameraParams[5];
	
	double mdW;             // distortion model coeff
	Vector2f mvCenter;     // Pixel projection center
	Vector2f mvFocal;      // Pixel focal length
	Vector2f mvImageSize;  
	float params[4];
	
	//resulting intermediary parameters
	Vector2f mvInvFocal;   // Inverse pixel focal length

	//for frustrum matrix computation
	Vector2f mvImplaneTL;   
	Vector2f mvImplaneBR;
	
};



