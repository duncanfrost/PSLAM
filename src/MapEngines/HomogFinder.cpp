#include "HomogFinder.h"
#include <Eigen/SVD>
#include <Eigen/LU>

float mdMaxPixelErrorSquared=25.;
float HomogScore(Matrix3f &Homog,p_match &_hmatch)
{
	double denom=Homog(2,0)*_hmatch.u1p+Homog(2,1)*_hmatch.v1p+Homog(2,2);
	Vector2f v2err;
	v2err[0]=(Homog(0,0)*_hmatch.u1p+Homog(0,1)*_hmatch.v1p+Homog(0,2))/denom -_hmatch.u1c;
	v2err[1]=(Homog(1,0)*_hmatch.u1p+Homog(1,1)*_hmatch.v1p+Homog(1,2))/denom -_hmatch.v1c;
	
	Vector2f v2errPix=v2err;
	float residue= v2errPix.transpose()*v2errPix;
	return residue;
}

float HomogScoreTukey(Matrix3f &Homog,p_match &_hmatch)
{
	float residue= HomogScore(Homog,_hmatch);
	
	if(residue<mdMaxPixelErrorSquared)
		return residue;
	else 
		return mdMaxPixelErrorSquared;
}
void HomographyFromMatchesRobustLM(std::vector<p_match> &_hmatches,Matrix3f &Homography)
{
	if(_hmatches.size()>4)
	{
		float mLMLambda = 0.0001;//before LevMarConstant
		float mdLambdaFactor = 2.0;

		for(int iter=0;iter<10;iter++)
		{

			//get sigma Tuckey
			std::vector<float> vdErrorsEssential;
			for(unsigned int i=0; i<_hmatches.size(); i++)
			{
				Vector2f x=Vector2f(_hmatches[i].u1p,_hmatches[i].v1p);
				Vector2f x_warped=homogWarp(x,Homography);
				Vector2f x_warped_des=Vector2f(_hmatches[i].u1c,_hmatches[i].v1c);
				float error=(x_warped_des-x_warped).transpose()*(x_warped_des-x_warped);
				vdErrorsEssential.push_back(error);
				//error+=(x_warped_des-x_warped).transpose()*(x_warped_des-x_warped);
			}
			float sigmaTukey=getSigmaSquared(vdErrorsEssential);
			if(sigmaTukey<100)sigmaTukey=100;
			
			float residue=0;
			MatrixXf H = MatrixXf::Zero(9, 9);//hessian
			MatrixXf G = MatrixXf::Zero(9, 1);//Jacobian*error
			
			for(unsigned int i=0; i<_hmatches.size(); i++)
			{
				Vector2f x=Vector2f(_hmatches[i].u1p,_hmatches[i].v1p);
				Vector2f x_warped=homogWarp(x,Homography);
				Vector2f x_warped_des=Vector2f(_hmatches[i].u1c,_hmatches[i].v1c);
				float error=(x_warped_des-x_warped).transpose()*(x_warped_des-x_warped);
				
				float TukeyCoef=squareRootTukey(error,sigmaTukey);
				if(TukeyCoef>0)
				{
					MatrixXf xwJac=JacHomogWarp(x,x_warped,Homography);
					
					H+=TukeyCoef*xwJac.transpose()*xwJac;
					G+=TukeyCoef*xwJac.transpose()*(x_warped_des-x_warped);
					residue+=TukeyCoef*error;
				}
			}
						
			//LM approach
			for(int i=0;i<9;i++)
				H(i,i)=(1.+mLMLambda)*H(i,i);
			
			Eigen::FullPivLU<MatrixXf> lu(H);
			
			//get vector representing Homography
			VectorXf p(9);
			p=getVectorFromHomog(Homography);
			
			p+=1.*(lu.inverse()*G);
			
			Matrix3f newHomography=getHomogFromVector(p);
			
			//check if any good
			float residueAfter=0;
			
			for(unsigned int i=0; i<_hmatches.size(); i++)
			{
				Vector2f x=Vector2f(_hmatches[i].u1p,_hmatches[i].v1p);
				Vector2f x_warped=homogWarp(x,Homography);
				Vector2f x_warped_des=Vector2f(_hmatches[i].u1c,_hmatches[i].v1c);
				float error=(x_warped_des-x_warped).transpose()*(x_warped_des-x_warped);
				
				float TukeyCoef=squareRootTukey(error,sigmaTukey);
				if(TukeyCoef>0)
					residueAfter=TukeyCoef*error;
			}
			
			if(residueAfter<residue)
			{
			  	mdLambdaFactor = 2.0;
				mLMLambda *= 0.3;
				Homography=newHomography;
			}
			else
			{
				mLMLambda = mLMLambda * mdLambdaFactor;
				mdLambdaFactor = mdLambdaFactor * 2;
			}
		}	
	}
}

Matrix3f HomographyFromMatchesRANSAC(std::vector<p_match> &_hmatches,int nb_trial)
{
	Matrix3f mm3BestHomography;
	if(_hmatches.size() < 10)
	{
		mm3BestHomography = HomographyFromMatches(_hmatches);
		return mm3BestHomography;
	}
	
	// Enough matches? Run MLESAC.
	int anIndices[4];

	mm3BestHomography = Matrix3f::Identity();
	double dBestError = 999999999999999999.9;

	// Do 500 MLESAC trials.
	for(int nR = 0; nR < nb_trial ; nR++)
	{ 
		// Find set of four unique matches
		for(int i=0; i<4; i++)
		{
			bool isUnique = false;
			int n;
			while(!isUnique)
			{
				n = rand() % _hmatches.size();
				isUnique =true;
				for(int j=0; j<i && isUnique; j++)
				if(anIndices[j] == n)
					isUnique = false;
			};
			anIndices[i] = n;
		}
		std::vector<p_match> vMinimalMatches;
		for(int i=0; i<4; i++)
			vMinimalMatches.push_back(_hmatches[anIndices[i]]);

		// Find a homography from the minimal set..
		Matrix3f m3Homography = HomographyFromMatches(vMinimalMatches);

		//..and sum resulting MLESAC score
		double dError = 0.0;
		for(unsigned int i=0; i<_hmatches.size(); i++)
			dError += HomogScoreTukey(m3Homography, _hmatches[i]);

		if(dError < dBestError)
		{
			mm3BestHomography = m3Homography;
			dBestError = dError;
		}
	}
	return mm3BestHomography;
	
}
Matrix3f HomographyFromMatches(std::vector<p_match> &_hmatches)
{
	unsigned int nPoints = _hmatches.size();
	if(nPoints>=4)
	{
		int nRows = 2*nPoints;
		if(nRows < 9)
			nRows = 9;
		MatrixXf m2Nx9(nRows,9);
		
		for(unsigned int n=0; n<nPoints; n++)
		{
			double u = _hmatches[n].u1c;
			double v = _hmatches[n].v1c;
			
			double x = _hmatches[n].u1p;
			double y = _hmatches[n].v1p;
			
			// [u v]T = H [x y]T
			m2Nx9(n*2+0,0) = x;
			m2Nx9(n*2+0,1) = y;
			m2Nx9(n*2+0,2) = 1;
			m2Nx9(n*2+0,3) = 0;
			m2Nx9(n*2+0,4) = 0;
			m2Nx9(n*2+0,5) = 0;
			m2Nx9(n*2+0,6) = -x*u;
			m2Nx9(n*2+0,7) = -y*u;
			m2Nx9(n*2+0,8) = -u;

			m2Nx9(n*2+1,0) = 0;
			m2Nx9(n*2+1,1) = 0;
			m2Nx9(n*2+1,2) = 0;
			m2Nx9(n*2+1,3) = x;
			m2Nx9(n*2+1,4) = y;
			m2Nx9(n*2+1,5) = 1;
			m2Nx9(n*2+1,6) = -x*v;
			m2Nx9(n*2+1,7) = -y*v;
			m2Nx9(n*2+1,8) = -v;
		}
		if(nRows == 9)  
		for(int i=0; i<9; i++)  // Zero the last row of the matrix, 
			m2Nx9(8,i) = 0.0;  // TooN SVD leaves out the null-space otherwise

		//JacobiSVD<MatrixXf> svdHomography(m2Nx9,ComputeFullV);
		JacobiSVD<MatrixXf> svdHomography(m2Nx9,ComputeThinV);
		const Eigen::MatrixXf V = svdHomography.matrixV();
		
		Matrix3f m3Homography;
		for(int i=0; i<3; i++) 
			for(int j=0; j<3; j++) 
				m3Homography(i,j)=V(3*i+j,8);
		
		return m3Homography;
	}
	else
		return Matrix3f::Identity();
}



static double SampsonusError(Vector2f &v2Dash, const Matrix3f &m3Essential, Vector2f &v2,Camera *_ptCam)
{
  Vector3f v3Dash = _ptCam->UnProjectZ1(v2Dash);
  Vector3f v3 = _ptCam->UnProjectZ1(v2);  
  
  double dError = v3Dash.dot(m3Essential * v3);
  
  Vector3f fv3 = m3Essential * v3;
  Vector3f fTv3Dash = m3Essential.transpose() * v3Dash;
  
  Vector2f fv3Slice; fv3Slice[0]=fv3[0];fv3Slice[1]=fv3[1];
  Vector2f fTv3DashSlice; fTv3DashSlice[0]= fTv3Dash[0];fTv3DashSlice[1]= fTv3Dash[1];
  
  return (dError * dError / (fv3Slice.dot(fv3Slice) + fTv3DashSlice.dot(fTv3DashSlice)));
}


VectorXf getVectorFromHomog(Matrix3f _homog)
{
	VectorXf res(9);
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			if(i==j)
				res[i+3*j]=_homog(i,j)-1;
			else
				res[i+3*j]=_homog(i,j);			
		}
	return res;
}
Matrix3f getHomogFromVector(VectorXf _v)
{
	Matrix3f res;
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			if(i==j)
				res(i,j)=_v[i+3*j]+1;
			else
				res(i,j)=_v[i+3*j];
		}
	return res;	
}

Vector2f homogWarp(Vector2f &x,Matrix3f &Homog)
{
	Vector2f res;
	double denom=1./(Homog(2,0)*x[0]+Homog(2,1)*x[1]+Homog(2,2));
	res[0]=(Homog(0,0)*x[0]+Homog(0,1)*x[1]+Homog(0,2))*denom;
	res[1]=(Homog(1,0)*x[0]+Homog(1,1)*x[1]+Homog(1,2))*denom;
	return res;
}
MatrixXf JacHomogWarp(Vector2f &x,Vector2f &xw,Matrix3f &Homog)
{
	MatrixXf dW(2,9);dW.setZero();
	double Denom=1./(Homog(2,0)*x[0]+Homog(2,1)*x[1]+Homog(2,2));
	dW(0,0)=x[0]*Denom;dW(0,2)=-x[0]*xw[0]*Denom;dW(0,3)=x[1]*Denom;dW(0,5)=-x[1]*xw[0]*Denom;dW(0,6)=Denom;dW(0,8)=-xw[0]*Denom;
	dW(1,1)=x[0]*Denom;dW(1,2)=-x[0]*xw[1]*Denom;dW(1,4)=x[1]*Denom;dW(1,5)=-x[1]*xw[1]*Denom;dW(1,7)=Denom;dW(1,8)=-xw[1]*Denom;
	return dW;
}
