#include "GlobalTransfoEstim.h"

cv::Mat warpImageInv(cv::Mat &Img,Matrix3f _homography)
{
	//cv::InputArray M;
	cv::Matx33f M(_homography(0,0),_homography(0,1),_homography(0,2),
		      _homography(1,0),_homography(1,1),_homography(1,2),
		      _homography(2,0),_homography(2,1),_homography(2,2)
	);
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	cv::warpPerspective(Img, res, M, res.size(), cv::INTER_LINEAR + cv::WARP_INVERSE_MAP);
  	//Eigen::FullPivLU<MatrixXf> lu(_homography);		
	//Matrix3f inv_homography=lu.inverse();
	/*
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	for(int x=0;x<Img.size().width;x++)
		for(int y=0;y<Img.size().height;y++)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			
			double denom=_homography(2,0)*x+_homography(2,1)*y+_homography(2,2);
			Vector2f pos2;
			pos2[0]=(_homography(0,0)*x+_homography(0,1)*y+_homography(0,2))/denom;
			pos2[1]=(_homography(1,0)*x+_homography(1,1)*y+_homography(1,2))/denom;
			
			//std::cout<<x<<","<<y<<" = > "<<pos2.transpose()<<std::endl;
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}*/
	return res;			
}
float getOverlapFromHomography(cv::Size imageSize,Matrix3f _homography,short modulo)
{
	int nb_in=0;
	int nb_sample=0;
	for(int x=0;x<imageSize.width;x+=modulo)
		for(int y=0;y<imageSize.height;y+=modulo)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			
			double denom=_homography(2,0)*x+_homography(2,1)*y+_homography(2,2);
			Vector2f pos2;
			pos2[0]=(_homography(0,0)*x+_homography(0,1)*y+_homography(0,2))/denom;
			pos2[1]=(_homography(1,0)*x+_homography(1,1)*y+_homography(1,2))/denom;
			
			if(pos2[0]>=0 && pos2[0]<=imageSize.width-1 && pos2[1]>=0 && pos2[1]<=imageSize.height-1)
				nb_in++;
			nb_sample++;
		}
	return (float)nb_in/nb_sample;	
}




Vector3f estimateRotationPyr(cv::Mat &T0,cv::Mat &I0,Camera* _cam,int pyr_lvl,int modulo)
{
	cv::Mat Ttemp[pyr_lvl];
	cv::Mat Itemp[pyr_lvl];
	cv::Mat *pT=&T0;
	cv::Mat *pI=&I0;
	

	if(pyr_lvl==1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );	
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
	else if(pyr_lvl>1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );					
		
		for(int l=1;l<pyr_lvl;l++)
		{
			cv::pyrDown( Ttemp[l-1], Ttemp[l], cv::Size( Ttemp[l-1].cols/2, Ttemp[l-1].rows/2 ) );
			cv::pyrDown( Itemp[l-1], Itemp[l], cv::Size( Itemp[l-1].cols/2, Itemp[l-1].rows/2 ) );					
		}
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
	
	
	cv::Mat &T=*pT;
	cv::Mat &I=*pI;

	//blur images
	/*int KernelSize=3;
	cv::GaussianBlur( T, T, cv::Size( KernelSize, KernelSize ), 0, 0 );
	cv::GaussianBlur( I, I, cv::Size( KernelSize, KernelSize ), 0, 0 );*/
	//cv::blur( T, T, cv::Size( KernelSize, KernelSize ), cv::Point(-1,-1) );
	//cv::blur( I, I, cv::Size( KernelSize, KernelSize ), cv::Point(-1,-1) );

	//initial rotation estimation
	Vector3f w; w.setZero();

	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		
		Matrix3f RotMat=getRotationFromThetaU(w);
		
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(ScaleLevel(pyr_lvl)* posInTemplate);
			Vector2f posInImage=applyRotationPyr(posInTemplate,w,_cam,pyr_lvl);
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width && posInImage[1]>=0 && posInImage[1]<I.size().height)
			{			
			
				
				Vector2f gradI;	
				gradI[0]=invScaleLevel(pyr_lvl)*getColorSubpixf(ImgGradx,posInImage);
				gradI[1]=invScaleLevel(pyr_lvl)*getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
				
				/*Vector2i posInImageE=toVector2i(posInImage);
				Vector2f gradI;	
				gradI[0]=invScaleLevel(pyr_lvl)*ImgGradx.at<float>(posInImageE[1],posInImageE[0]);
				gradI[1]=invScaleLevel(pyr_lvl)*ImgGrady.at<float>(posInImageE[1],posInImageE[0]);
				
				//color in current image
				float colorI=I.at<uchar>(posInImageE[1],posInImageE[0]);*/
				
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
				MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
				MatrixXf JIw(1,3);JIw=gradI.transpose()*dProj_dpos2*x2_cross;
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
		
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt(Dw.squaredNorm())<0.0005 || residue==0)
		{
			//std::cout<<"nb iter = "<<iter<<std::endl;
			break;
		}
	}
	return w;
}

Vector3f estimateRotation(cv::Mat &T,cv::Mat &I,Camera* _cam,int modulo)
{
	//initial rotation estimation
	Vector3f w; w.setZero();

	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		
		Matrix3f RotMat=getRotationFromThetaU(w);
		
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(posInTemplate);
			Vector2f posInImage=applyRotation(posInTemplate,w,_cam);
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width && posInImage[1]>=0 && posInImage[1]<I.size().height)
			{			
			
				Vector2f gradI;	gradI[0]=getColorSubpixf(ImgGradx,posInImage);gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
				
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
				MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
				MatrixXf JIw(1,3);JIw=gradI.transpose()*dProj_dpos2*x2_cross;
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
		
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt(Dw.squaredNorm())<0.0001 || residue==0)
			break;
	}
	return w;
}

Vector3f estimateRotation(std::vector<p_match> &matches, Camera* _cam)
{
	//std::cout<<"estimateRotation"<<std::endl;
	//initial rotation estimation
	Vector3f w; w.setZero();


	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		
		//just to monitor evolution
		float residue=0;
		
		
		Matrix3f RotMat=getRotationFromThetaU(w);
		
		for(int i=0;i<matches.size();i++)
		{
			Vector2f posInTemplate;	posInTemplate[0]=matches[i].u1p;posInTemplate[1]=matches[i].v1p;
			
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(posInTemplate);
			Vector2f posInImage=applyRotation(posInTemplate,w,_cam);
			
			Vector2f posDes=Vector2f(matches[i].u1c,matches[i].v1c);
			
			Vector2f error=posInImage-posDes;
			
			MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
			MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
			MatrixXf Jw(2,3);Jw=dProj_dpos2*x2_cross;
			
			//update derivatives
			Jte+=Jw.transpose()*error;
			H+=Jw.transpose()*Jw;
			
			residue+=sqrt(error.transpose()*error);

		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
		
		
		//std::cout<<iter<<" : normalized residue = "<<residue/matches.size()<<std::endl;
		if(sqrt(Dw.squaredNorm())<0.0001 || residue==0)
			break;
	}
	return w;
}
Vector3f estimateRotationInv(std::vector<p_match> &matches, Camera* _cam)
{
	//std::cout<<"estimateRotation"<<std::endl;
	//initial rotation estimation
	Vector3f w; w.setZero();


	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		
		//just to monitor evolution
		//float residue=0;
		
		
		Matrix3f RotMat=getRotationFromThetaU(w);
		
		for(int i=0;i<matches.size();i++)
		{
			Vector2f posInTemplate;	posInTemplate[0]=matches[i].u1c;posInTemplate[1]=matches[i].v1c;
			
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(posInTemplate);
			Vector2f posInImage=_cam->Project(posMeter2);

			
			Vector2f posDes=Vector2f(matches[i].u1p,matches[i].v1p);
			
			Vector2f error=posInImage-posDes;
			
			MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
			MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
			MatrixXf Jw(2,3);Jw=dProj_dpos2*x2_cross;
			
			//update derivatives
			Jte+=Jw.transpose()*error;
			H+=Jw.transpose()*Jw;
			
			//residue+=sqrt(error.transpose()*error);

		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
		
		
		//std::cout<<iter<<" : normalized residue = "<<residue/matches.size()<<std::endl;
		if(sqrt(Dw.squaredNorm())<0.0001)
			break;
	}
	return w;
}
Vector3f SE2toSO3(Vector3f _p,Camera* _cam)
{
	//transform two points with SE2 and do non linear optim to find same thing using SO3
	Vector2f refPoints[2];
	refPoints[0]=_cam->getCenterImage()/2;
	refPoints[1]=3.*_cam->getCenterImage()/2;
	
	//transform using SE
	Vector2f trefPoints_d[2];
	for(int i=0;i<2;i++)
		trefPoints_d[i]=applySE2(refPoints[i],_p,_cam);
	
	//do non linear optim;
	Vector3f w;	w.setZero();
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		Matrix3f RotMat=getRotationFromThetaU(w);
		for(int i=0;i<2;i++)
		{
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(refPoints[i]);
			Vector2f trefPoints=applyRotation(refPoints[i],w,_cam);
			
			Vector2f error = trefPoints-trefPoints_d[i];
			
			residue+=error.transpose()*error;
			
			MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
			MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
			MatrixXf JIw(2,3);JIw=dProj_dpos2*x2_cross;
			
			//update derivatives
			Jte+=JIw.transpose()*error;
			H+=JIw.transpose()*JIw;
			
		}
		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
				
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(residue<0.1)
			break;
	}
	return w;
}

Vector3f estimateSE2(cv::Mat &T,cv::Mat &I,Camera* _cam,int modulo)
{
	//initial rotation estimation
	Vector3f p; p.setZero();

	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			Vector2f posInImage=applySE2(posInTemplate,p,_cam);
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width && posInImage[1]>=0 && posInImage[1]<I.size().height)
			{			
			
				Vector2f gradI;	gradI[0]=getColorSubpixf(ImgGradx,posInImage);gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
				
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf JIw(1,3);
				JIw(0,0)=gradI[0];JIw(1,0)=gradI[1];
				JIw(2,0)=gradI.transpose()*getSO2RotationFromAplha(p[2]+M_PI/2.)*(posInTemplate-_cam->getCenterImage());
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dp=-lambda*Jte;
		
		if(!isnan(Dp[0]))
			p+=Dp;
		
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt( (Dp.segment(0,2).squaredNorm())<1. && Dp[2]<0.01)|| residue==0)
			break;
	}
	return p;
}

#define AMOVERBOSE
Vector3f estimateSE2Pyr(cv::Mat &T0,cv::Mat &I0,Camera* _cam,int pyr_lvl,int modulo,int _max_iter)
{
	cv::Mat Ttemp[pyr_lvl];
	cv::Mat Itemp[pyr_lvl];
	cv::Mat *pT=&T0;
	cv::Mat *pI=&I0;
	

	if(pyr_lvl==1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );	
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
	else if(pyr_lvl>1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );					
		
		for(int l=1;l<pyr_lvl;l++)
		{
			cv::pyrDown( Ttemp[l-1], Ttemp[l], cv::Size( Ttemp[l-1].cols/2, Ttemp[l-1].rows/2 ) );
			cv::pyrDown( Itemp[l-1], Itemp[l], cv::Size( Itemp[l-1].cols/2, Itemp[l-1].rows/2 ) );					
		}
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
		
	cv::Mat &T=*pT;
	cv::Mat &I=*pI;

	//blur images
	/*int KernelSize=3;
	cv::GaussianBlur( T, T, cv::Size( KernelSize, KernelSize ), 0, 0 );
	cv::GaussianBlur( I, I, cv::Size( KernelSize, KernelSize ), 0, 0 );*/
	/*int KernelSize=3;
	cv::blur( T, T, cv::Size( KernelSize, KernelSize ), cv::Point(-1,-1) );
	cv::blur( I, I, cv::Size( KernelSize, KernelSize ), cv::Point(-1,-1) );*/

	//initial rotation estimation
	Vector3f w; w.setZero();
	
	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	
	//iterative non linear optimisation
	for(int iter=0;iter<_max_iter;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		
		//#pragma omp parallel for shared(Jte, H, residue) private(pyr_lvl)
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			Vector2f posInImage=applySE2Pyr(posInTemplate,w,_cam,pyr_lvl);
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width-1 && posInImage[1]>=0 && posInImage[1]<I.size().height-1)
			{			
			
				
				Vector2f gradI;	
				gradI[0]=getColorSubpixf(ImgGradx,posInImage);
				gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
								
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf JIw(1,3);
				JIw(0,0)=invScaleLevel(pyr_lvl)*gradI[0];
				JIw(0,1)=invScaleLevel(pyr_lvl)*gradI[1];
				JIw(0,2)=invScaleLevel(pyr_lvl)*
				gradI.transpose()*getSO2RotationFromAplha(w[2]+M_PI/2.)*(ScaleLevel(pyr_lvl)*posInTemplate-_cam->getCenterImage());
				
//invScaleLevel(_lvl)*(getSO2RotationFromAplha(p[2])*(ScaleLevel(_lvl)*pos-_cam->getCenterImage())+_cam->getCenterImage()+p.segment(0,2));

				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}


		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=-lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
		
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt( (Dw.segment(0,2).squaredNorm())<1. && Dw[2]<0.01)|| residue==0)
		{
			//std::cout<<"nb iter = "<<iter<<std::endl;
			break;
		}
	}
	return w;
}

cv::Mat translateImage(cv::Mat &Img,Vector2f t)
{
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	for(int x=0;x<Img.size().width;x++)
		for(int y=0;y<Img.size().height;y++)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			Vector2f pos2;pos2=pos-t;
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}
	return res;
			
}

Vector2f estimateTranslation(cv::Mat &T,cv::Mat &I, int modulo)
{
	//initial rotation estimation
	Vector2f t; t.setZero();

	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);
	
	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(2);Jte.setZero();
		MatrixXf H(2,2);H.setZero();
		float residue=0;
		
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			Vector2f posInImage;posInImage=posInTemplate+t;
			if(posInImage[0]>=0 && posInImage[0]<I.size().width && posInImage[1]>=0 && posInImage[1]<I.size().height)
			{			
				Vector2f gradI;	gradI[0]=getColorSubpixf(ImgGradx,posInImage);gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
				
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=(float)colorT-colorI;
				residue+=error*error;
				
				MatrixXf JIw(1,2);JIw=gradI.transpose();
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}
		
		
		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector2f Dt=lambda*Jte;
		
		t+=Dt;
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt(Dt.squaredNorm())<0.05 || residue==0)
			break;

	}
	return t;
}
Vector2f applyRotation(Vector2f pos,Vector3f w,Camera* _cam)
{
	Vector3f posMeter2=getRotationFromThetaU(w)*_cam->UnProjectZ1(pos);
	return _cam->Project(posMeter2);
}
Matrix2f getSO2RotationFromAplha(float _a)
{
	Matrix2f res;
	res(0,0)=cos(_a);res(0,1)=-sin(_a);
	res(1,0)=sin(_a);res(1,1)=cos(_a);
	return res;
}
Vector3f invSE2(Vector3f p)
{
	Vector3f res;
	res.segment(0,2)=-getSO2RotationFromAplha(-p[2])*p.segment(0,2);
	res[2]=-p[2];
	return res;
}
Vector2f applySE2(Vector2f pos,Vector3f p,Camera* _cam)
{
	return getSO2RotationFromAplha(p[2])*(pos-_cam->getCenterImage())+_cam->getCenterImage()+p.segment(0,2);
}
Vector2f applySE2Pyr(Vector2f pos,Vector3f p,Camera* _cam,int _lvl)
{
	return invScaleLevel(_lvl)*(getSO2RotationFromAplha(p[2])*(ScaleLevel(_lvl)*pos-_cam->getCenterImage())+_cam->getCenterImage()+p.segment(0,2));
}

Vector2f applyRotationPyr(Vector2f pos,Vector3f w,Camera* _cam,int _lvl)
{
	Vector3f posMeter2=getRotationFromThetaU(w)*_cam->UnProjectZ1(ScaleLevel(_lvl)*  pos);
	return invScaleLevel(_lvl)*_cam->Project(posMeter2);
}
cv::Mat rotateImage(cv::Mat &Img,Vector3f w,Camera* _cam)
{
	
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	for(int x=0;x<Img.size().width;x++)
		for(int y=0;y<Img.size().height;y++)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			Vector2f pos2=applyRotation(pos,-w,_cam);
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}
	return res;			
}

cv::Mat rotateImageKeepAll(cv::Mat &Img,Vector3f w,Camera* _cam,Vector2f &offsetRes)
{
	//warp for corners of image
	Vector2f Corners1[4];
	Corners1[0][0]=0;Corners1[0][1]=0;
	Corners1[1][0]=Img.size().width;Corners1[1][1]=0;
	Corners1[2][0]=Img.size().width;Corners1[2][1]=Img.size().height;
	Corners1[3][0]=0;Corners1[3][1]=Img.size().height;
  
	Vector2f CornersRot[4];
	for(int i=0;i<4;i++)
		CornersRot[i]=applyRotation(Corners1[i],w,_cam);
	
	//get supporting x,y axis
	float minx,maxx;
	float miny,maxy;
	
	minx=CornersRot[0][0];maxx=CornersRot[0][0];
	miny=CornersRot[0][1];maxy=CornersRot[0][1];
	
	for(int i=1;i<4;i++)
	{
		if(minx>CornersRot[i][0])minx=CornersRot[i][0];
		if(miny>CornersRot[i][1])miny=CornersRot[i][1];
		if(maxx<CornersRot[i][0])maxx=CornersRot[i][0];
		if(maxy<CornersRot[i][1])maxy=CornersRot[i][1];
	}

	//will place our new 0,0 to minx,minx
	offsetRes[0]=(int)(minx-0.5);
	offsetRes[1]=(int)(miny-0.5);
	
	int new_width=(int)(maxx+0.5)-offsetRes[0];
	int new_height=(int)(maxy+0.5)-offsetRes[1];
  
	/*cv::Mat res;res.create(new_height, new_width, CV_8UC1);
	for(int x=0;x<new_width;x++)
		for(int y=0;y<new_height;y++)
		{
			Vector2f pos;pos[0]=x+offsetRes[0];pos[1]=y+offsetRes[1];
			Vector2f pos2=applyRotation(pos,-w,_cam);
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}*/
		
	//a bit faster
	Matrix3f RotationMat=getRotationFromThetaU(-w);
	cv::Mat res;res.create(new_height, new_width, CV_8UC1);
	for(int x=0;x<new_width;x++)
		for(int y=0;y<new_height;y++)
		{
			Vector2f pos;pos[0]=x+offsetRes[0];pos[1]=y+offsetRes[1];
			Vector3f posMeter2=RotationMat*_cam->UnProjectZ1(pos);
			Vector2f pos2=_cam->Project(posMeter2);
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}
		
	return res;			
}
cv::Mat applySE2Image(cv::Mat &Img,Vector3f p,Camera* _cam)
{
	
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	Vector3f invp=invSE2(p);
	for(int x=0;x<Img.size().width;x++)
		for(int y=0;y<Img.size().height;y++)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			Vector2f pos2=applySE2(pos,invp,_cam);
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}
	return res;			
}

Vector3f estimateRotationInverse(cv::Mat &T,cv::Mat &I,Camera* _cam,int modulo)
{
	//initial rotation estimation
	Vector3f w; w.setZero();

	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(T,ImgGradx,ImgGrady);

	//iterative non linear optimisation
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		
		Matrix3f RotMat=getRotationFromThetaU(w);
		
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(posInTemplate);
			Vector2f posInImage=applyRotation(posInTemplate,w,_cam);
			if(posInImage[0]>=0 && posInImage[0]<I.size().width && posInImage[1]>=0 && posInImage[1]<I.size().height)
			{			
			
				Vector2f gradI;	gradI[0]=getColorSubpixf(ImgGradx,posInImage);gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(T,posInImage);
				
				//color in template:
				uchar colorT=I.at<uchar>(y,x);
				
				float error=(float)colorT-colorI;
				residue+=error*error;
				
				MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
				MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
				MatrixXf JIw(1,3);JIw=gradI.transpose()*dProj_dpos2*x2_cross;
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
			}
		}

		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		w+=-Dw;
		
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(sqrt(Dw.squaredNorm())<0.0005 || residue==0)
			break;
	}
	return -w;
}


Vector4f invSE2Tz(Vector4f p)
{
	Vector4f res;
	res.segment(0,2)=-(1./(1.+p[3]))*getSO2RotationFromAplha(-p[2])*p.segment(0,2);
	res[2]=-p[2];
	res[3]=1./(1.+p[3])-1.;
	return res;
}
Vector2f applySE2Tz(Vector2f pos,Vector4f p,Camera* _cam)
{
	return (1+p[3])*getSO2RotationFromAplha(p[2])*(pos-_cam->getCenterImage())+_cam->getCenterImage()+p.segment(0,2);
}

Vector2f applySE2TzPyr(Vector2f pos,Vector4f p,Camera* _cam,int _lvl)
{
	return invScaleLevel(_lvl)*((1+p[3])*getSO2RotationFromAplha(p[2])*(ScaleLevel(_lvl)*pos-_cam->getCenterImage())+_cam->getCenterImage()+p.segment(0,2));
}

cv::Mat applySE2TzImage(cv::Mat &Img,Vector4f p,Camera* _cam)
{
	
	cv::Mat res;res.create(Img.size().height, Img.size().width, CV_8UC1);
	Vector4f invp=invSE2Tz(p);
	for(int x=0;x<Img.size().width;x++)
		for(int y=0;y<Img.size().height;y++)
		{
			Vector2f pos;pos[0]=x;pos[1]=y;
			Vector2f pos2=applySE2Tz(pos,invp,_cam);
			
			if(pos2[0]>=0 && pos2[0]<=Img.size().width-1 && pos2[1]>=0 && pos2[1]<=Img.size().height-1)
				res.at<uchar>(y,x)=getColorSubpix(Img,pos2);
			else
				res.at<uchar>(y,x)=0;
		}
	return res;			
}

HomogeneousMatrix SE2TztoSO3(Vector4f _p,float mean_depth,Camera* _cam)
{
	//transform two points with SE2 and do non linear optim to find same thing using SO3
	Vector2f refPoints[2];
	refPoints[0]=_cam->getCenterImage()/2;
	refPoints[1]=3.*_cam->getCenterImage()/2;
	
	//transform using SE
	/*Vector2f trefPoints_d[2];
	for(int i=0;i<2;i++)
		trefPoints_d[i]=applySE2Tz(refPoints[i],_p,_cam);
	
	//do non linear optim;
	Vector4f p;	p.setZero();
	Vector3f w;float t_z;
	for(int iter=0;iter<50;iter++)
	{
		w=p.segment(0,3);
		t_z=p[3];
		
		VectorXf Jte(4);Jte.setZero();
		MatrixXf H(4,4);H.setZero();
		float residue=0;
		Matrix3f RotMat=getRotationFromThetaU(w);
		for(int i=0;i<2;i++)
		{
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(refPoints[i]);
			Vector2f trefPoints=applyRotation(refPoints[i],w,_cam);
			
			Vector2f error = trefPoints-trefPoints_d[i];
			
			residue+=error.transpose()*error;
			
			MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
			MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
			MatrixXf JIw(2,4);
			JIw.block(0,0,2,3)=dProj_dpos2*x2_cross;
			//JIw.block(0,3,2,1)=...; TODO
			
			//update derivatives
			Jte+=JIw.transpose()*error;
			H+=JIw.transpose()*JIw;
			
		}
		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dp=lambda*Jte;
		
		if(!isnan(Dp[0]))
			p+=Dp;
				
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(residue<0.02)
			break;
	}*/
	Vector2f trefPoints_d[2];
	for(int i=0;i<2;i++)
		trefPoints_d[i]=applySE2(refPoints[i],_p.segment(0,3),_cam);
	
	//do non linear optim;
	Vector3f w;	w.setZero();
	for(int iter=0;iter<50;iter++)
	{
		VectorXf Jte(3);Jte.setZero();
		MatrixXf H(3,3);H.setZero();
		float residue=0;
		Matrix3f RotMat=getRotationFromThetaU(w);
		for(int i=0;i<2;i++)
		{
			Vector3f posMeter2=RotMat*_cam->UnProjectZ1(refPoints[i]);
			Vector2f trefPoints=applyRotation(refPoints[i],w,_cam);
			
			Vector2f error = trefPoints-trefPoints_d[i];
			
			residue+=error.transpose()*error;
			
			MatrixXf dProj_dpos2(2,3);dProj_dpos2=_cam->m2PixProjJac()*_cam->ProjectZ1_Jac_X(posMeter2);
			MatrixXf x2_cross(3,3);x2_cross=GetSkew(posMeter2);
			MatrixXf JIw(2,3);JIw=dProj_dpos2*x2_cross;
			
			//update derivatives
			Jte+=JIw.transpose()*error;
			H+=JIw.transpose()*JIw;
			
		}
		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector3f Dw=lambda*Jte;
		
		if(!isnan(Dw[0]))
			w+=Dw;
				
		//std::cout<<iter<<" : residue = "<<residue<<std::endl;
		if(residue<0.1)
			break;
	}
	
	return HomogeneousMatrix(0,0,mean_depth/(1+_p[3])-mean_depth,w[0],w[1],w[2]);
}

float estimateSE2TzPyr(cv::Mat &T0,cv::Mat &I0,Camera* _cam,Vector4f &res,int pyr_lvl,int modulo,int _max_iter)
{
	cv::Mat Ttemp[pyr_lvl];
	cv::Mat Itemp[pyr_lvl];
	cv::Mat *pT=&T0;
	cv::Mat *pI=&I0;
	

	if(pyr_lvl==1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );	
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
	else if(pyr_lvl>1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );					
		
		for(int l=1;l<pyr_lvl;l++)
		{
			cv::pyrDown( Ttemp[l-1], Ttemp[l], cv::Size( Ttemp[l-1].cols/2, Ttemp[l-1].rows/2 ) );
			cv::pyrDown( Itemp[l-1], Itemp[l], cv::Size( Itemp[l-1].cols/2, Itemp[l-1].rows/2 ) );					
		}
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
		
	cv::Mat &T=*pT;
	cv::Mat &I=*pI;

	//initial rotation estimation
	Vector4f p; p.setZero();
	
	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	
	//iterative non linear optimisation
	float residue=0;
	int nb_valid_pix=0;
	for(int iter=0;iter<_max_iter;iter++)
	{
		VectorXf Jte(4);Jte.setZero();
		MatrixXf H(4,4);H.setZero();
		residue=0;
		nb_valid_pix=0;
		
		//#pragma omp parallel for shared(Jte, H, residue) private(pyr_lvl)
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			Vector2f posInImage=applySE2TzPyr(posInTemplate,p,_cam,pyr_lvl);
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width-1 && posInImage[1]>=0 && posInImage[1]<I.size().height-1)
			{			
			
				
				Vector2f gradI;	
				gradI[0]=getColorSubpixf(ImgGradx,posInImage);
				gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
								
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf JIw(1,4);
				JIw(0,0)=invScaleLevel(pyr_lvl)*gradI[0];
				JIw(0,1)=invScaleLevel(pyr_lvl)*gradI[1];
				JIw(0,2)=invScaleLevel(pyr_lvl)*
				gradI.transpose()*getSO2RotationFromAplha(p[2]+M_PI/2.)*(ScaleLevel(pyr_lvl)*posInTemplate-_cam->getCenterImage());
			
				JIw(0,3)=invScaleLevel(pyr_lvl)*gradI.transpose()*getSO2RotationFromAplha(p[2])*(ScaleLevel(pyr_lvl)*posInTemplate-_cam->getCenterImage());
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
				nb_valid_pix++;
			}
		}


		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector4f Dp=-lambda*Jte;
		
		if(!isnan(Dp[0]))
			p+=Dp;
		
		
		//std::cout<<iter<<" : residue = "<<residue/nb_valid_pix<<std::endl;
		if(sqrt( (Dp.segment(0,2).squaredNorm())<1. && Dp[2]<0.01 && Dp[3]<0.005)|| residue==0)
		{
			//std::cout<<"nb iter = "<<iter<<std::endl;
			break;
		}
	}
	res=p;
	//return residue/nb_valid_pix;
	
	//return ZMSSD
	int ISum = 0;
	int ISumSq = 0;
	int TSum = 0;
	int TSumSq = 0;
	int CrossSum = 0;
	int nb_valid=0;
	int nb_tempted=0;
	for(int x=0;x<T.size().width;x+=modulo)
		for(int y=0;y<T.size().height;y+=modulo)
	{
		Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
		Vector2f posInImage=applySE2TzPyr(posInTemplate,p,_cam,pyr_lvl);
		nb_tempted++;
		if(posInImage[0]>=0 && posInImage[0]<I.size().width-1 && posInImage[1]>=0 && posInImage[1]<I.size().height-1)
		{			

			//float I=getColorSubpix(I,posInImage);
			int It=I.at<uchar>((int)(posInImage[1]+0.5),(int)(posInImage[0]+0.5));
			int Tt=T.at<uchar>(y,x);
			ISum += It;
			ISumSq += It*It;
			TSum += Tt;
			TSumSq += Tt*Tt;
			CrossSum += It * Tt;
			nb_valid++;
		}

	}
	if(nb_valid<0.7*nb_tempted || p[3]>0.3 || p[3]<-0.3)//if too much out or too much tz estimated then considered as diverged
	{
		//if(nb_valid<0.7*I.size().width*I.size().height)std::cout<<"ratio vis = "<<(float)nb_valid/nb_tempted<<std::endl;
		//else if (p[3]>0.3 || p[3]<-0.3)std::cout<<"tz = "<<p[3]<<std::endl;
		
		return 1e10;
	}
	else	
		return ((2*TSum*ISum - TSum*TSum - ISum*ISum)/nb_valid + ISumSq + TSumSq - 2*CrossSum)/nb_valid;
	
}
float estimateTzPyr(cv::Mat &T0,cv::Mat &I0,Camera* _cam,int pyr_lvl,int modulo,int _max_iter)
{
	cv::Mat Ttemp[pyr_lvl];
	cv::Mat Itemp[pyr_lvl];
	cv::Mat *pT=&T0;
	cv::Mat *pI=&I0;
	

	if(pyr_lvl==1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );	
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
	else if(pyr_lvl>1)
	{
		cv::pyrDown( T0, Ttemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );
		cv::pyrDown( I0, Itemp[0], cv::Size( T0.size().width/2,T0.size().height/2 ) );					
		
		for(int l=1;l<pyr_lvl;l++)
		{
			cv::pyrDown( Ttemp[l-1], Ttemp[l], cv::Size( Ttemp[l-1].cols/2, Ttemp[l-1].rows/2 ) );
			cv::pyrDown( Itemp[l-1], Itemp[l], cv::Size( Itemp[l-1].cols/2, Itemp[l-1].rows/2 ) );					
		}
		pT=&Ttemp[pyr_lvl-1];
		pI=&Itemp[pyr_lvl-1];
	}
		
	cv::Mat &T=*pT;
	cv::Mat &I=*pI;

	VectorXf p(1);p[0]=0;
	
	//get current iamge gradients
	cv::Mat ImgGradx,ImgGrady;
	imageGradients(I,ImgGradx,ImgGrady);

	
	//iterative non linear optimisation
	float residue=0;
	int nb_valid_pix=0;
	for(int iter=0;iter<_max_iter;iter++)
	{
		VectorXf Jte(1);Jte.setZero();
		MatrixXf H(1,1);H.setZero();
		residue=0;
		nb_valid_pix=0;
		
		//#pragma omp parallel for shared(Jte, H, residue) private(pyr_lvl)
		for(int x=0;x<T.size().width;x+=modulo)
			for(int y=0;y<T.size().height;y+=modulo)
		{
			Vector2f posInTemplate;	posInTemplate[0]=x;posInTemplate[1]=y;
			Vector2f posInImage=invScaleLevel(pyr_lvl)*((1+p[0])*(ScaleLevel(pyr_lvl)*posInTemplate-_cam->getCenterImage())+_cam->getCenterImage());
			
			if(posInImage[0]>=0 && posInImage[0]<I.size().width-1 && posInImage[1]>=0 && posInImage[1]<I.size().height-1)
			{			
			
				
				Vector2f gradI;	
				gradI[0]=getColorSubpixf(ImgGradx,posInImage);
				gradI[1]=getColorSubpixf(ImgGrady,posInImage);
				
				//color in current image
				float colorI=getColorSubpix(I,posInImage);
								
				//color in template:
				uchar colorT=T.at<uchar>(y,x);
				
				float error=colorI-(float)colorT;
				residue+=error*error;
				
				MatrixXf JIw(1,1);
				JIw(0,0)=invScaleLevel(pyr_lvl)*gradI.transpose()*getSO2RotationFromAplha(p[2])*(ScaleLevel(pyr_lvl)*posInTemplate-_cam->getCenterImage());
				
				//update derivatives
				Jte+=JIw.transpose()*error;
				H+=JIw.transpose()*JIw;
				nb_valid_pix++;
			}
		}


		float lambda=1.;
		H.ldlt().solveInPlace(Jte);
		Vector4f Dp=-lambda*Jte;
		
		if(!isnan(Dp[0]))
			p+=Dp;
		
		
		//std::cout<<iter<<" : residue = "<<residue/nb_valid_pix<<std::endl;

	}
	return p[0];
}