//Homogemeous matrix implementation


#pragma once

#include <objectr3d/AmoDefines.h>
#include <objectr3d/GeoFunctions.h>
#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;

#define PISDF M_PI
#define PIdiv180SDF (PISDF/180.0)

class HomogeneousMatrix
{
public:
	HomogeneousMatrix();
	//init from minimal representation rodriguez style : t first and rotation u
	HomogeneousMatrix(float *_p);
	HomogeneousMatrix(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz);
	HomogeneousMatrix(VectorXf _p);
	//init from translation and rotation matrix
	HomogeneousMatrix(Matrix4f _Hmat);
	HomogeneousMatrix(Vector3f _tVec,Matrix3f _Rmat);
	
	void Init(float *_p);
	void Init(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz);
	void Init(VectorXf _p);
	void Init(Matrix4f _Hmat);
	void Init(Vector3f _tVec,Matrix3f _Rmat);
	void Init(HomogeneousMatrix &_H);
	void SetIdentity();
	
	
	//get full representation
	Vector3f get_translation(){return t;};
	Matrix3f get_rotation(){return Rmat;}
	Matrix4f get_HomogMatrix(){return Hmat;};
	
	//get 6D parameterisation
	void get_p(float *_p){for(int i=0;i<3;i++)_p[i]=u[i];for(int i=0;i<3;i++)_p[i+3]=w[i];}
	VectorXf get_p(){VectorXf p(6);for(int i=0;i<3;i++)p[i]=u[i];for(int i=0;i<3;i++)p[i+3]=w[i];return p;}
	Vector3f get_u(){return u;};
	Vector3f get_w(){return w;};
	float get_angle(){return angle;};

	//get inverse of homogeneous transformation
	HomogeneousMatrix inverse();
	
	//set functions
	void set_translation(Vector3f _tVec){Init(_tVec,Rmat);};
	void set_rotation(Matrix3f _Rmat){Init(t,_Rmat);};
	
	//rescale translational part
	void normalize_translation(float _d){set_translation((_d/sqrt(t.squaredNorm()))*t);};

	Vector3f operator*(const Vector3f &_p);
	HomogeneousMatrix operator*(HomogeneousMatrix _H);
	HomogeneousMatrix operator+(HomogeneousMatrix _H);
	HomogeneousMatrix operator/(float div);
	
	friend std::ostream& operator<<(std::ostream& os, HomogeneousMatrix dt);
	void printMatrix44();
	
	//function to apply simple movement to camera pose
	void RotateX ( float Angle );
	void RotateY ( float Angle );
	void RotateZ ( float Angle );
	void TranslateX ( float Distance );
	void TranslateY ( float Distance );
	void TranslateZ ( float Distance );
	
	//input output file functions
	void saveToFile(char *filename);
	bool loadFromFile(char *filename);
	void saveToStream(std::ofstream &fout);
	void loadFromStream(std::ifstream &fout);
	
	~HomogeneousMatrix(){};	
private:
	//homogeneous representation
	Matrix4f Hmat;
	Matrix3f Rmat;
	Vector3f t;
	
	//minimal representation
	Vector3f u;
	Vector3f w;	
	float angle;
	
};



