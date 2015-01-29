#include "HomogeneousMatrix.h"

HomogeneousMatrix22::HomogeneousMatrix22()
{
    std::cout << sizeof(*this) << std::endl;
    std::cout << sizeof(HomogeneousMatrix22) << std::endl;
    Matrix4f _Hmat=Matrix4f::Identity();
	Init(_Hmat);
}
HomogeneousMatrix22::HomogeneousMatrix22(float *_p)
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=_p[i];
	Init(p);
}
void HomogeneousMatrix22::SetIdentity()
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=0;
	Init(p);
}

HomogeneousMatrix22::HomogeneousMatrix22(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz)
{
	VectorXf p(6);
	p[0]=_tx;p[1]=_ty;p[2]=_tz;
	p[3]=_wx;p[4]=_wy;p[5]=_wz;
	Init(p);
}
HomogeneousMatrix22::HomogeneousMatrix22(VectorXf _p)
{
	Init(_p);
}
HomogeneousMatrix22::HomogeneousMatrix22(Matrix4f _Hmat)
{
	Init(_Hmat);
}
HomogeneousMatrix22::HomogeneousMatrix22(Vector3f _tVec,Matrix3f _Rmat)
{
	Init(_tVec,_Rmat);
}

void HomogeneousMatrix22::Init(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz)
{
	VectorXf p(6);
	p[0]=_tx;p[1]=_ty;p[2]=_tz;
	p[3]=_wx;p[4]=_wy;p[5]=_wz;
	Init(p);
}

void HomogeneousMatrix22::Init(float *_p)
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=_p[i];
	Init(p);
}
void HomogeneousMatrix22::Init(VectorXf _p)
{

}
void HomogeneousMatrix22::Init(Matrix4f _Hmat)
{

	
}
void HomogeneousMatrix22::Init(Vector3f _tVec,Matrix3f _Rmat)
{

	
}


void HomogeneousMatrix22::Init(HomogeneousMatrix22 &_H)
{
	Init(_H.get_HomogMatrix());
}
	

Vector3f HomogeneousMatrix22::operator*(const Vector3f &_p)
{
    Matrix3f Rmat;
    Vector3f t;
	return Rmat*_p+t;
}
HomogeneousMatrix22 HomogeneousMatrix22::operator*(HomogeneousMatrix22 _H2)
{
    Matrix3f Rmat;
        Vector3f t;
    HomogeneousMatrix22 res(Rmat*_H2.get_translation()+t,Rmat*_H2.get_rotation());
	return res;
}
HomogeneousMatrix22 HomogeneousMatrix22::operator+(HomogeneousMatrix22 _H2)
{
	VectorXf p(6);p=_H2.get_p()+get_p();
    HomogeneousMatrix22 res(p);
	return res;
}
HomogeneousMatrix22 HomogeneousMatrix22::operator/(float div)
{
	VectorXf p(6);p=get_p()/div;
    HomogeneousMatrix22 res(p);
	return res;
}
/*HomogeneousMatrix HomogeneousMatrix::operator+(HomogeneousMatrix _H2)
{
	HomogeneousMatrix res(Rmat*_H2.get_translation()+t,Rmat*_H2.get_rotation(),true);
	return res;
}*/
HomogeneousMatrix22 HomogeneousMatrix22::inverse()
{
	//Matrix3f Rinv=getRotationFromThetaU(-1.*w);
    Matrix3f Rmat;
	Eigen::FullPivLU<MatrixXf> lu(Rmat);
	Matrix3f Rinv =lu.inverse();	
	
        Vector3f t;

    HomogeneousMatrix22 res(-1.*Rinv*t,Rinv);
	return res;
}

std::ostream& operator<<(std::ostream& os, HomogeneousMatrix22 dt)
{
    os << dt.get_u().transpose()<<", "<< dt.get_w().transpose();
    return os;
}
void HomogeneousMatrix22::printMatrix44()
{
    Matrix4f Hmat;
	std::cout<<Hmat<<std::endl;
}
void HomogeneousMatrix22::RotateX ( float Angle )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
    velocity.setZero();velocity[3]=-Angle*0.017;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix22::RotateY ( float Angle )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
    velocity.setZero();velocity[4]=-Angle*0.017;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix22::RotateZ ( float Angle )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
    velocity.setZero();velocity[5]=-Angle*0.017;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};

void HomogeneousMatrix22::TranslateX ( float Distance )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
	velocity.setZero();velocity[0]=-Distance;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};

void HomogeneousMatrix22::TranslateY ( float Distance )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
	velocity.setZero();velocity[1]=-Distance;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix22::TranslateZ ( float Distance )
{
     Matrix4f Hmat;
	VectorXf velocity(6);
	velocity.setZero();velocity[2]=-Distance;
    HomogeneousMatrix22 c2_to_c(velocity);
    HomogeneousMatrix22 poseCurrent(Hmat);
    HomogeneousMatrix22 newH=c2_to_c*poseCurrent;
	Init(newH);
};

#include <fstream>
void HomogeneousMatrix22::saveToFile(char *filename)
{

}
bool HomogeneousMatrix22::loadFromFile(char *filename)
{
	std::ifstream fout;
	fout.open(filename);
	if(!fout.is_open())
	{
		std::cerr<<"HomogeneousMatrix : problem loading file "<<filename<<std::endl;
		return false;
	}
	VectorXf p(6);
	for(int j=0;j<6;j++)
		fout.read((char*)&p[j],sizeof(float));
	Init(p);
	fout.close();	
	return true;
}

/*void HomogeneousMatrix::saveToStream(std::ofstream &fout)
{
	for(int j=0;j<3;j++)fout.write((const char*)&u[j],sizeof(float));
	for(int j=0;j<3;j++)fout.write((const char*)&w[j],sizeof(float));
}
void HomogeneousMatrix::loadFromStream(std::ifstream &fout)
{
	VectorXf p(6);
	for(int j=0;j<6;j++)fout.read((char*)&p[j],sizeof(float));
	Init(p);
}*/
/*void HomogeneousMatrix::saveToStream(std::ofstream &fout)
{
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
		fout.write((const char*)&Hmat(i,j),sizeof(float));
}
void HomogeneousMatrix::loadFromStream(std::ifstream &fout)
{
	Matrix4f lHmat;
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			fout.read((char*)&lHmat(i,j),sizeof(float));
	Init(lHmat);
}*/
void HomogeneousMatrix22::saveToStream(std::ofstream &fout)
{


	
}
void HomogeneousMatrix22::loadFromStream(std::ifstream &fout)
{

}
