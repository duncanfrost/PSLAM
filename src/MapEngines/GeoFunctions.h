#pragma once

#include "../AmoDefines.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <Eigen/Geometry>

using namespace Eigen;

Matrix3f GetSkew(Vector3f _w);
Vector3f GetInvSkew(Matrix3f _mat);
Matrix4f GetHomogeneousMatrix(float *_p);
void LnSe3(Matrix4f _H,float *_p);
Matrix3f getRotationFromThetaU(Vector3f _w);
Vector3f getThetaUFromRotation(Matrix3f _R);

Vector3f toHomogeneous(Vector2f _x);
Vector4f toHomogeneous(Vector3f _x);

bool inImage(Vector2f _p,int _width,int _height);

/*struct PointMatch
{
	PointMatch(){id=-1;};
	
	Vector2f v2CamPlaneFirst;
	Vector2f v2CamPlaneSecond;
	int id;
	int idGT;
};*/


float gaussianNoise();
//Vector2i toVector2i(Vector2f _v);
//Vector2f toVector2f(Vector2i _v);

#include <vector>

//TODO will have to put that in MapOptimiser.h
//non null elements of jacobian matrix corresponding to the derivative of 
//reprojection error with respect to point pos dx and camera pose dp
struct kf_vis_jacobian
{
	Vector2f proj_error;//projection error
	int pt_index;//index of pt in optimisation list
	MatrixXf de_dx;//jacobian wrt point position
	int cam_index;//index of cam in optimisation list
	MatrixXf de_dp;//jacobian wrt cam position
	float weight;
};

//#define MINSIGMATUKEY 0.16 //minimum value in BA
#define MINSIGMATUKEY 2. //minimum value in BA
float getSigmaSquared(std::vector<float> &vdErrorSquared);
float squareRootTukey(float errorSquared, float sigma);

Matrix3f differentiateNormalisedVectorByItself(Vector3f v);

