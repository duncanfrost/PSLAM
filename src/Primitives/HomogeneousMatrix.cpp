#include "HomogeneousMatrix.h"

HomogeneousMatrix::HomogeneousMatrix()
{
	Matrix4f _Hmat=Matrix4f::Identity();
	Init(_Hmat);
}
HomogeneousMatrix::HomogeneousMatrix(float *_p)
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=_p[i];
	Init(p);
}
void HomogeneousMatrix::SetIdentity()
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=0;
	Init(p);
}

HomogeneousMatrix::HomogeneousMatrix(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz)
{
	VectorXf p(6);
	p[0]=_tx;p[1]=_ty;p[2]=_tz;
	p[3]=_wx;p[4]=_wy;p[5]=_wz;
	Init(p);
}
HomogeneousMatrix::HomogeneousMatrix(VectorXf _p)
{
	Init(_p);
}
HomogeneousMatrix::HomogeneousMatrix(Matrix4f _Hmat)
{
	Init(_Hmat);
}
HomogeneousMatrix::HomogeneousMatrix(Vector3f _tVec,Matrix3f _Rmat)
{
	Init(_tVec,_Rmat);
}

void HomogeneousMatrix::Init(float _tx,float _ty,float _tz,float _wx,float _wy,float _wz)
{
	VectorXf p(6);
	p[0]=_tx;p[1]=_ty;p[2]=_tz;
	p[3]=_wx;p[4]=_wy;p[5]=_wz;
	Init(p);
}

void HomogeneousMatrix::Init(float *_p)
{
	VectorXf p(6);
	for(int i=0;i<6;i++)p[i]=_p[i];
	Init(p);
}
void HomogeneousMatrix::Init(VectorXf _p)
{
	for(int i=0;i<3;i++)u[i]=_p[i];
	for(int i=0;i<3;i++)w[i]=_p[i+3];
	angle=sqrt(w.squaredNorm());
	
	Hmat.setZero();	
		
	if(angle!=0)
	{
		AngleAxis<float> AngAx(angle,w/angle);
		Rmat=AngAx.toRotationMatrix();
	}
	else
		Rmat=Matrix3f::Identity();
	

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			Hmat(i,j)=Rmat(i,j);

	//set translation
	Matrix3f SkewW=GetSkew(w);
	Matrix3f A=Matrix3f::Identity()+SkewW*(1-cos(angle))/(angle*angle)+SkewW*SkewW*(angle-sin(angle))/(angle*angle*angle);
	if(isnan(A(0,0)))A=Matrix3f::Identity();
	
	t=A*u;
	for(int i=0;i<3;i++)Hmat(i,3)=t[i];

	Hmat(3,3)=1;
}
void HomogeneousMatrix::Init(Matrix4f _Hmat)
{
	Hmat=_Hmat;
	Hmat.block(0,0,3,1)=Hmat.block(0,0,3,1)/sqrt(Hmat.block(0,0,3,1).squaredNorm());
	Hmat.block(0,1,3,1)=Hmat.block(0,1,3,1)/sqrt(Hmat.block(0,1,3,1).squaredNorm());
	Hmat.block(0,2,3,1)=Hmat.block(0,2,3,1)/sqrt(Hmat.block(0,2,3,1).squaredNorm());
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			Rmat(i,j)=Hmat(i,j);
	for(int i=0;i<3;i++)t[i]=Hmat(i,3);
	
	AngleAxis<float> AngAx;
	AngAx=AngAx.fromRotationMatrix(Rmat);
	angle=AngAx.angle();
	w=AngAx.axis()*angle;
	
	Matrix3f skewW=GetSkew(w);
	Matrix3f Ainv=Matrix3f::Identity()-0.5*skewW+skewW*skewW*(2*sin(angle)-angle*(1+cos(angle)))/(2.*angle*angle*sin(angle));
	if(isnan(Ainv(0,0)))Ainv=Matrix3f::Identity();
	
	u=Ainv*t;
	
}
void HomogeneousMatrix::Init(Vector3f _tVec,Matrix3f _Rmat)
{
	Rmat=_Rmat;
	Rmat.block(0,0,3,1)=Rmat.block(0,0,3,1)/sqrt(Rmat.block(0,0,3,1).squaredNorm());
	Rmat.block(0,1,3,1)=Rmat.block(0,1,3,1)/sqrt(Rmat.block(0,1,3,1).squaredNorm());
	Rmat.block(0,2,3,1)=Rmat.block(0,2,3,1)/sqrt(Rmat.block(0,2,3,1).squaredNorm());

	t=_tVec;
	
	Hmat.setZero();
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			Hmat(i,j)=Rmat(i,j);
	for(int i=0;i<3;i++)Hmat(i,3)=t[i];
	Hmat(3,3)=1;
	
	AngleAxis<float> AngAx;
	AngAx=AngAx.fromRotationMatrix(Rmat);
	angle=AngAx.angle();
	w=AngAx.axis()*angle;
	//std::cout<<"HomogeneousMatrix >> angle = "<<angle<<std::endl;
	
	Matrix3f skewW=GetSkew(w);
	Matrix3f Ainv=Matrix3f::Identity()-0.5*skewW+skewW*skewW*(2*sin(angle)-angle*(1+cos(angle)))/(2.*angle*angle*sin(angle));
	if(isnan(Ainv(0,0)))Ainv=Matrix3f::Identity();


	u=Ainv*t;
	
}


void HomogeneousMatrix::Init(HomogeneousMatrix &_H)
{
	Init(_H.get_HomogMatrix());
}
	

Vector3f HomogeneousMatrix::operator*(const Vector3f &_p)
{
	return Rmat*_p+t;
}
HomogeneousMatrix HomogeneousMatrix::operator*(HomogeneousMatrix _H2)
{
	HomogeneousMatrix res(Rmat*_H2.get_translation()+t,Rmat*_H2.get_rotation());
	return res;
}
HomogeneousMatrix HomogeneousMatrix::operator+(HomogeneousMatrix _H2)
{
	VectorXf p(6);p=_H2.get_p()+get_p();
	HomogeneousMatrix res(p);
	return res;
}
HomogeneousMatrix HomogeneousMatrix::operator/(float div)
{
	VectorXf p(6);p=get_p()/div;
	HomogeneousMatrix res(p);
	return res;
}
/*HomogeneousMatrix HomogeneousMatrix::operator+(HomogeneousMatrix _H2)
{
	HomogeneousMatrix res(Rmat*_H2.get_translation()+t,Rmat*_H2.get_rotation(),true);
	return res;
}*/
HomogeneousMatrix HomogeneousMatrix::inverse()
{
	//Matrix3f Rinv=getRotationFromThetaU(-1.*w);
	
	Eigen::FullPivLU<MatrixXf> lu(Rmat);
	Matrix3f Rinv =lu.inverse();	
	
	HomogeneousMatrix res(-1.*Rinv*t,Rinv);
	return res;
}

std::ostream& operator<<(std::ostream& os, HomogeneousMatrix dt)
{
    os << dt.get_u().transpose()<<", "<< dt.get_w().transpose();
    return os;
}
void HomogeneousMatrix::printMatrix44()
{
	std::cout<<Hmat<<std::endl;
}
void HomogeneousMatrix::RotateX ( float Angle )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[3]=-Angle*PIdiv180SDF;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix::RotateY ( float Angle )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[4]=-Angle*PIdiv180SDF;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix::RotateZ ( float Angle )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[5]=-Angle*PIdiv180SDF;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};

void HomogeneousMatrix::TranslateX ( float Distance )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[0]=-Distance;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};

void HomogeneousMatrix::TranslateY ( float Distance )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[1]=-Distance;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};
void HomogeneousMatrix::TranslateZ ( float Distance )
{
	VectorXf velocity(6);
	velocity.setZero();velocity[2]=-Distance;
	HomogeneousMatrix c2_to_c(velocity);
	HomogeneousMatrix poseCurrent(Hmat);
	HomogeneousMatrix newH=c2_to_c*poseCurrent;
	Init(newH);
};

#include <fstream>
void HomogeneousMatrix::saveToFile(char *filename)
{
	std::ofstream fout;
	fout.open(filename);
	
	for(int j=0;j<3;j++)
		fout.write((const char*)&u[j],sizeof(float));
	for(int j=0;j<3;j++)
		fout.write((const char*)&w[j],sizeof(float));
	fout.close();
}
bool HomogeneousMatrix::loadFromFile(char *filename)
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
void HomogeneousMatrix::saveToStream(std::ofstream &fout)
{
	char isPTAMversion=0;
	fout.write((const char*)&isPTAMversion,sizeof(char));
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
		fout.write((const char*)&Hmat(i,j),sizeof(float));
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		fout.write((const char*)&Rmat(i,j),sizeof(float));
	for(int i=0;i<3;i++)
		fout.write((const char*)&t[i],sizeof(float));
	for(int i=0;i<3;i++)
		fout.write((const char*)&u[i],sizeof(float));
	for(int i=0;i<3;i++)
		fout.write((const char*)&w[i],sizeof(float));
	fout.write((const char*)&angle,sizeof(float));
	
}
void HomogeneousMatrix::loadFromStream(std::ifstream &fout)
{
	char isPTAMversion;
	fout.read((char*)&isPTAMversion,sizeof(char));//in PTAM when save transformation, only save H matrix
	if(isPTAMversion==0)
	{
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				fout.read((char*)&Hmat(i,j),sizeof(float));
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				fout.read((char*)&Rmat(i,j),sizeof(float));
		for(int i=0;i<3;i++)
			fout.read((char*)&t[i],sizeof(float));
		for(int i=0;i<3;i++)
			fout.read((char*)&u[i],sizeof(float));
		for(int i=0;i<3;i++)
			fout.read((char*)&w[i],sizeof(float));
		fout.read((char*)&angle,sizeof(float));
	}
	else
	{
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				fout.read((char*)&Hmat(i,j),sizeof(float));
			
		Init(Hmat);
	}
}
