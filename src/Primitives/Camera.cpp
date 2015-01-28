#include "Camera.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>


Camera::Camera(int _width, int _height, CameraType _camType)
{
	Init(_width,_height,_camType);
}
Camera Camera::getHalfSizeCamera()
{
	Camera resCam;
	float *params1=getParametersPTAM();
	float *params2=resCam.getParametersPTAM();
	
	//need to keep same vector mgvvCameraParams
	for(int i=0;i<5;i++)params2[i]=params1[i];
	
	//Init cam with half size=> all parameters are set properly
	resCam.InitParams(mvImageSize[0]/2,mvImageSize[1]/2);
	
	return resCam;
}

void Camera::Init(int _width, int _height, CameraType _camType)
{
 	bool isCalibrated = false;
	std::string calib_filename;
	switch (_camType)
	{
		case CamPlaystationEye:
		{
			calib_filename=PEYE_CALIB_FILE;
			
			isCalibrated=true;
		}
		break;
		
		case KittiCam:
		{
			calib_filename=KITTI_CALIB_FILE;
			
			isCalibrated=true;
		}
		break;
		
		case CamFlea2:
		{
			calib_filename="/media/Data/Datas/Calibrations/PointGrey/camera.cfg";
			isCalibrated=true;
		}
		break;

		default:
		break;
	}
	
	if(!isCalibrated)
	{
		mgvvCameraParams[0]=0.5;mgvvCameraParams[1]= 0.75;mgvvCameraParams[2]= 0.5;mgvvCameraParams[3]= 0.5;mgvvCameraParams[4]= 0.1;
	}
	else
	{
		//read calib_filename and get params
		std::ifstream in1;	
		in1.open(calib_filename.c_str());
		if (in1.is_open())
		{
			int taille_buf=500;		char buffer[taille_buf];		in1.getline(buffer, taille_buf); 
			std::string mot;mot=buffer;for(int i=0;i<20;i++)mot[i]=' ';
			std::istringstream ligne;		ligne.str(mot);
			ligne>>mgvvCameraParams[0]>>mgvvCameraParams[1]>>mgvvCameraParams[2]>>mgvvCameraParams[3]>>mgvvCameraParams[4];;
			//std::cout<<"mgvvCameraParams = "<<mgvvCameraParams[0]<<" "<<mgvvCameraParams[1]<<" "<<mgvvCameraParams[2]<<" "<<mgvvCameraParams[3]<<" "<<mgvvCameraParams[4]<<std::endl;	
			in1.close();
		}
		else
		{
			coutRed<<"Calib file not found"<<endlRed;
			mgvvCameraParams[0]=0.5;mgvvCameraParams[1]= 0.75;mgvvCameraParams[2]= 0.5;mgvvCameraParams[3]= 0.5;mgvvCameraParams[4]= 0.1;
		}
	}
	InitParams(_width,_height);
}
		
void Camera::InitParams(int _width, int _height)
{
	
	mvImageSize[0]=_width;
	mvImageSize[1]=_height;
	
	// First: Focal length and image center in pixel coordinates
	mvFocal[0] = mvImageSize[0] * mgvvCameraParams[0];
	mvFocal[1] = mvImageSize[1] * mgvvCameraParams[1];
	mvCenter[0] = mvImageSize[0] * mgvvCameraParams[2] - 0.5;
	mvCenter[1] = mvImageSize[1] * mgvvCameraParams[3] - 0.5;
	
	params[0]=mvCenter[0];
	params[1]=mvCenter[1];
	params[2]=mvFocal[0];
	params[3]=mvFocal[1];

	// Some radial distortion parameters..
	mdW =  mgvvCameraParams[4];	

	// One over focal length
	mvInvFocal[0] = 1.0 / mvFocal[0];
	mvInvFocal[1] = 1.0 / mvFocal[1];
  
	// Work out the linear projection values for the UFB
	{
		// First: Find out how big the linear bounding rectangle must be
		std::vector<Vector2f > vv2Verts;
		vv2Verts.push_back(UnProject(Vector2f( -0.5, -0.5)));
		vv2Verts.push_back(UnProject(Vector2f( mvImageSize[0]-0.5, -0.5)));
		vv2Verts.push_back(UnProject(Vector2f( mvImageSize[0]-0.5, mvImageSize[1]-0.5)));
		vv2Verts.push_back(UnProject(Vector2f( -0.5, mvImageSize[1]-0.5)));
		Vector2f v2Min = vv2Verts[0];
		Vector2f v2Max = vv2Verts[0];
		for(int i=0; i<4; i++)
		for(int j=0; j<2; j++)
		{
			if(vv2Verts[i][j] < v2Min[j]) v2Min[j] = vv2Verts[i][j];
			if(vv2Verts[i][j] > v2Max[j]) v2Max[j] = vv2Verts[i][j];
		}
		mvImplaneTL = v2Min;
		mvImplaneBR = v2Max;
	}
}

Vector2f Camera::Project(const Vector2f& mvLastDistCam)
{
  Vector2f v2Im;
  v2Im[0] = mvLastDistCam[0]*mvFocal[0] + mvCenter[0];
  v2Im[1] = mvLastDistCam[1]*mvFocal[1] + mvCenter[1];
  return v2Im;
}

Vector2f Camera::ToMeters(const Vector2f& v2Im)
{
  Vector2f mvLastDistCam;
  mvLastDistCam[0] =(v2Im[0]-mvCenter[0])/mvFocal[0];
  mvLastDistCam[1] =(v2Im[1]-mvCenter[1])/mvFocal[1];
  return mvLastDistCam;	
}

Vector2f Camera::Project(const Vector3f& mvLastDistCam)
{
  Vector2f v2Im;
  v2Im[0] = mvLastDistCam[0]/mvLastDistCam[2]*mvFocal[0] + mvCenter[0];
  v2Im[1] = mvLastDistCam[1]/mvLastDistCam[2]*mvFocal[1] + mvCenter[1];
  return v2Im;
}
Vector2f Camera::ProjectZ1(const Vector3f& mvLastDistCam)
{
  Vector2f v2Im;
  v2Im[0] = mvLastDistCam[0]/mvLastDistCam[2];
  v2Im[1] = mvLastDistCam[1]/mvLastDistCam[2];
  return v2Im;
}
MatrixXf Camera::ProjectZ1_Jac_X(const Vector3f& mvLastDistCam)
{
  MatrixXf res(2,3);
  float invZ=1./mvLastDistCam[2];
  res(0,0)= invZ;res(0,1)= 0.;res(0,2)= -mvLastDistCam[0]*invZ*invZ;
  res(1,0)= 0.;res(1,1)= invZ;res(1,2)= -mvLastDistCam[1]*invZ*invZ;
  return res;
}
MatrixXf Camera::ProjectZ1_Jac_Dp(const Vector3f& mvLastDistCam)
{
  MatrixXf res(2,6);
  float invZ=1./mvLastDistCam[2];
  float x=mvLastDistCam[0]*invZ;
  float y=mvLastDistCam[1]*invZ;
  //float x=mvLastDistCam[0];
  //float y=mvLastDistCam[1];
  res(0,0)= invZ;res(0,1)= 0.;res(0,2)= -x*invZ;res(0,3)= -x*y;res(0,4)= (1+x*x);res(0,5)= -y;
  res(1,0)= 0.;res(1,1)= invZ;res(1,2)= -y*invZ;res(1,3)= -(1+y*y);res(1,4)= x*y;res(1,5)= x;
  return res;
}
//unproject linear by default at the moment 
//TODO : integrate distortion
Vector2f Camera::UnProject(const Vector2f& v2Im)
{
  Vector2f mvLastDistCam;
  mvLastDistCam[0] = (v2Im[0] - mvCenter[0]) * mvInvFocal[0];
  mvLastDistCam[1] = (v2Im[1] - mvCenter[1]) * mvInvFocal[1];
  return mvLastDistCam;
}
Vector3f Camera::UnProjectZ1(const Vector2f& v2Im)
{
  Vector3f mvLastDistCam;
  mvLastDistCam[0] = (v2Im[0] - mvCenter[0]) * mvInvFocal[0];
  mvLastDistCam[1] = (v2Im[1] - mvCenter[1]) * mvInvFocal[1];
  mvLastDistCam[2] = 1.;
  return mvLastDistCam;
}
Matrix2f Camera::m2PixProjJac()
{
	Matrix2f res;
	res(0,0)=mvFocal[0];res(0,1)=0;
	res(1,0)=0;res(1,1)=mvFocal[1];
	return res;	
}
Matrix2f Camera::Pix2mProjJac()
{
	Matrix2f res;
	res(0,0)=1./mvFocal[0];res(0,1)=0;
	res(1,0)=0;res(1,1)=1./mvFocal[1];
	return res;	
}
Matrix4f Camera::getGlProjectionMat(double near, double far)
{
	Matrix4f mproj;mproj.setZero();

	double left = mvImplaneTL[0] * near;
	double right = mvImplaneBR[0] * near;
	double top = mvImplaneTL[1] * near;
	double bottom = mvImplaneBR[1] * near;

	mproj(0,0) = (2 * near) / (right - left);
	mproj(1,1) = (2 * near) / (top - bottom);

	mproj(0,2) = (right + left) / (left - right);
	mproj(1,2) = (top + bottom) / (bottom - top);
	mproj(2,2) = (far + near) / (far - near);
	mproj(3,2) = 1;

	mproj(2,3) = 2*near*far / (near - far);

	return mproj;
}

