#include "GeoFunctions.h"



Matrix3f GetSkew(Vector3f _w)
{
	Matrix3f res;
	res.setZero();
			  res(0,1)=-_w[2];	res(0,2)=_w[1];
	res(1,0)=_w[2];				res(1,2)=-_w[0];
	res(2,0)=-_w[1];  res(2,1)=_w[0];
	return res;
}

Vector3f GetInvSkew(Matrix3f _mat)
{
	Vector3f res;
	res[0]=_mat(2,1);
	res[1]=_mat(0,2);
	res[2]=_mat(1,0);
	return res;
}
Matrix4f GetHomogeneousMatrix(float *_p)
{
	Vector3f u;for(int i=0;i<3;i++)u[i]=_p[i];
	Vector3f w;for(int i=0;i<3;i++)w[i]=_p[i+3];
	float angle=sqrt(w.squaredNorm());
	
	Matrix4f res;
	res.setZero();
	
		
	//set rotation
	Matrix3f rmat;
	//rmat=GetSkew(w).exp();
	Matrix3f SkewW=GetSkew(w);
	if(angle==0)	rmat=Matrix3f::Identity();
	else 		rmat=Matrix3f::Identity()+SkewW*sin(angle)/angle+SkewW*SkewW*(1-cos(angle))/(angle*angle);

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			res(i,j)=rmat(i,j);

	//set translation
	Matrix3f A;	
	if(angle==0)	A=Matrix3f::Identity();
	else 		A=Matrix3f::Identity()+SkewW*(1-cos(angle))/(angle*angle)+SkewW*SkewW*(angle-sin(angle))/(angle*angle*angle);
	Vector3f t=A*u;
	for(int i=0;i<3;i++)res(i,3)=t[i];

	res(3,3)=1;
	return res;
}

void LnSe3(Matrix4f _H,float *_p)
{
	Matrix3f rmat;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			rmat(i,j)=_H(i,j);
	Vector3f t;for(int i=0;i<3;i++)t[i]=_H(i,3);
	
		
	float phi=acos((rmat.trace()-1)/2.);
	Matrix3f skewW;
	if(phi==0)	skewW.setZero();
	else 		skewW=(rmat-rmat.transpose())*phi/(2.*sin(phi));
	Vector3f w=GetInvSkew(skewW);
	float angle=sqrt(w.squaredNorm());
	
	
	Matrix3f Ainv;
	
	if(angle==0)	Ainv=Matrix3f::Identity();
	else 		Ainv=Matrix3f::Identity()-0.5*skewW+skewW*skewW*(2*sin(angle)-angle*(1+cos(angle)))/(2.*angle*angle*sin(angle));
	Vector3f u;	u=Ainv*t;
	
	for(int i=0;i<3;i++)_p[i]=u[i];
	for(int i=0;i<3;i++)_p[i+3]=w[i];
	//std::cout<<rmat.log()<<std::endl;
	
}
Vector3f getThetaUFromRotation(Matrix3f Rot)
{
	AngleAxis<float> AngAx;
	AngAx=AngAx.fromRotationMatrix(Rot);
	float angle=AngAx.angle();
	return AngAx.axis()*angle;
}

bool inImage(Vector2f _p,int _width,int _height)
{
	if(_p[0]>0 && _p[1]>0 && _p[0]<_width-1 && _p[1]<_height-1)
		return true;
	else 
		return false;
}

Matrix3f getRotationFromThetaU(Vector3f w)
{
	/*float angle=sqrt(_w.squaredNorm());
	Matrix3f SkewW=GetSkew(_w);
	if(angle==0)	return Matrix3f::Identity();
	else 		return Matrix3f::Identity()+SkewW*sin(angle)/angle+SkewW*SkewW*(1-cos(angle))/(angle*angle);*/
	
	float angle=sqrt(w.squaredNorm());
	
	if(angle!=0)
	{
		AngleAxis<float> AngAx(angle,w/angle);
		return AngAx.toRotationMatrix();
	}
	else
		return Matrix3f::Identity();
}

Vector3f toHomogeneous(Vector2f _x)
{
	Vector3f res;
	res[0]=_x[0];res[1]=_x[1];res[2]=1.;
	return res;
}
Vector4f toHomogeneous(Vector3f _x)
{
	Vector4f res;
	res[0]=_x[0];res[1]=_x[1];res[2]=_x[2];res[3]=1.;	
	return res;
}


float gaussianNoise()
{
	//srand(time(0));
	double GaussNum = 0.0;
	int NumInSum = 10;
	for(int i = 0; i < NumInSum; i++)
	{
		GaussNum += ((double)rand()/(double)RAND_MAX - 0.5);
	}
	GaussNum = GaussNum*sqrt((double)12/(double)NumInSum);
	return GaussNum;
}
Vector2i toVector2i(Vector2f _v){return Vector2i(_v[0],_v[1]);};
Vector2f toVector2f(Vector2i _v){return Vector2f(_v[0],_v[1]);};


float getSigmaSquared(std::vector<float> &vdErrorSquared)
{ 
  float dSigmaSquared; 
  assert(vdErrorSquared.size() > 0);
  std::sort(vdErrorSquared.begin(), vdErrorSquared.end());
  float dMedianSquared = vdErrorSquared[vdErrorSquared.size() / 2];
  float dSigma = 1.4826 * (1 + 5.0 / (vdErrorSquared.size() * 2 - 6)) * sqrt(dMedianSquared);
  dSigma =  4.6851 * dSigma;
  dSigmaSquared = dSigma * dSigma;
  return dSigmaSquared;
}



float squareRootTukey(float errorSquared, float sigma)
{
	if(errorSquared > sigma)
		return 0.0;
	else
		return 1.0 - (errorSquared / sigma);
}
Matrix3f differentiateNormalisedVectorByItself(Vector3f v)
{
	float norm=sqrt(v.squaredNorm());
	Matrix3f res=(1./norm)*Matrix3f::Identity();
	res-=(1./(norm*norm*norm))*v*v.transpose();
	return res;
}
