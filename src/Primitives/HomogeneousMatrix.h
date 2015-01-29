//Homogemeous matrix implementation


#pragma once

#include "../AmoDefines.h"
#include "../MapEngines/GeoFunctions.h"
#include <Eigen/Core>
#include <iostream>

using namespace Eigen;


class HomogeneousMatrix22
{
public:
    HomogeneousMatrix22();
	//init from minimal representation rodriguez style : t first and rotation u
    HomogeneousMatrix22(float *_p);
    HomogeneousMatrix22(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz);
    HomogeneousMatrix22(VectorXf _p);
	//init from translation and rotation matrix
    HomogeneousMatrix22(Matrix4f _Hmat);
    HomogeneousMatrix22(Vector3f _tVec,Matrix3f _Rmat);
	
	void Init(float *_p);
	void Init(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz);
	void Init(VectorXf _p);
	void Init(Matrix4f _Hmat);
	void Init(Vector3f _tVec,Matrix3f _Rmat);
    void Init(HomogeneousMatrix22 &_H);
	void SetIdentity();
	
	
	//get full representation
    Vector3f get_translation(){Vector3f test; return test;};
    Matrix3f get_rotation(){Matrix3f Rmat; return Rmat;}
    Matrix4f get_HomogMatrix(){Matrix4f Hmat; return Hmat;};
	
	//get 6D parameterisation
    void get_p(float *_p){}
    VectorXf get_p(){VectorXf p(6); return p;}
    Vector3f get_u(){Vector3f u; return u;};
    Vector3f get_w(){Vector3f w; return w;};
    float get_angle(){float angle = 0; return angle;};

	//get inverse of homogeneous transformation
    HomogeneousMatrix22 inverse();
	
	//set functions
    void set_translation(Vector3f _tVec)
    {
        Matrix3f Rmat;
        Init(_tVec,Rmat);
    };
    void set_rotation(Matrix3f _Rmat){};
	
	//rescale translational part
    void normalize_translation(float _d){};

	Vector3f operator*(const Vector3f &_p);
    HomogeneousMatrix22 operator*(HomogeneousMatrix22 _H);
    HomogeneousMatrix22 operator+(HomogeneousMatrix22 _H);
    HomogeneousMatrix22 operator/(float div);
	
    friend std::ostream& operator<<(std::ostream& os, HomogeneousMatrix22 dt);
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
	
    ~HomogeneousMatrix22(){};
private:
	//homogeneous representation
    Matrix4f Hmat;
    float angle;




	
	//minimal representation


	
};



