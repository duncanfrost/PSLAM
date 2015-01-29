//perform BA using few small baseline images and no initial poses (random)
//this has been coded to test [3D Reconstruction from Accidental Motion - Fisher Yu]

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/MapEngines/MapOptimiser.h"
#include "../src/MapEngines/MapOptimiserEssential.h"
#include "../src/Primitives/obMap.h"

void Idle(void) ;

//Visualization
VisualisationModule *VisuEngine;

Camera myCamera;//camera object, calibration ...

//reprojection function and user interface
void keyboard_process(unsigned char key, int x, int y);
void addDrawFunction();
void disturbPose();
void disturbDepth();
void DoMiniBA();
void createRandomEstim();

#define nb_kf 5

EmptyWindow *InKfWindow;
MapWindow *mMapWindow;
PlottingWindow *PlotWindow;


HomogeneousMatrix PoseRef;
HomogeneousMatrix PoseEstim[nb_kf];
HomogeneousMatrix PoseGT[nb_kf];

float image_noise=0.;

int im_w;
int im_h;

//load pair of pose and 3D points from real scene
void createMap();
std::vector<Vector3f> RandomPoints;
std::vector<PointInvDepth> invDepthGT;
std::vector<PointInvDepth> invDepthEstim;

Vector2f getImagePlanNoise()
{
	Vector2f noise;
	noise[0]=gaussianNoise()*image_noise;
	noise[1]=gaussianNoise()*image_noise;	
	return myCamera.Pix2mProjJac()*noise;
}


int main(int argc, char** argv)
{
	im_w=640;
	im_h=480;
	myCamera.Init(im_w,im_h,CamPlaystationEye);
	
	//load GT map
	std::cout<<"load GT map"<<std::endl;
	createMap();
	createRandomEstim();
	//disturbPose();

	
	VisuEngine= new VisualisationModule(&Idle);
	//create a window to visualize reprojection and to grab keyboard event (see keyboard_process)
	VisuEngine->addWindowEmpty("Key frame 0 view",640,480,InKfWindow);
	VisuEngine->setOnDraw(addDrawFunction);
	VisuEngine->setOnKeyPress(keyboard_process);
	//create a window to visualize GT map and estimated one
	VisuEngine->addWindowMap("Map",640,480,&myCamera,mMapWindow);
	VisuEngine->setOnKeyPress(keyboard_process);
	
	VisuEngine->addWindowPlot("PoseError",640,480,PlotWindow,2,200,false);
	PlotWindow->setValuesPerPlot(0,1);//reproj error
	PlotWindow->setValuesPerPlot(1,2);//pose error
	
	//display gt map and estimated map
	mMapWindow->addCamera(&PoseRef,Vector3f(1,0,0),2);
	for(int i=0;i<nb_kf;i++)
	{
		mMapWindow->addCamera(&PoseGT[i],Vector3f(0,0,1),1);
		mMapWindow->addCamera(&PoseEstim[i],Vector3f(0,1,0),1);
	}
	
	mMapWindow->addPointCloud(&invDepthGT);
	mMapWindow->addPointCloud(&invDepthEstim);
	
	//HomogeneousMatrix moveCam(0,0.8,0.8,0.7,0,0);
	HomogeneousMatrix moveCam(0,2.1,2.1,2.0,0,0);
	mMapWindow->moveCamera(moveCam);
	
	
	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


void Idle(void) 
{
	float reproj_error=0;
 	PlotWindow->setVal(0,0,reproj_error);
 
	float error_trans=0;
	float error_rot=0;
	for(int i=0;i<nb_kf;i++)
	{
		HomogeneousMatrix errorPose=PoseGT[i]*PoseEstim[i].inverse();
		VectorXf error_p=errorPose.get_p();
		//std::cout<<"error_p = "<<error_p.transpose()<<std::endl;
		error_trans+=sqrt(error_p.segment(0,3).squaredNorm());
		error_rot+=sqrt(error_p.segment(3,3).squaredNorm());
	}

	PlotWindow->setVal(1,0,error_trans);
	PlotWindow->setVal(1,1,error_rot);
	
	PlotWindow->incrementTimeLine();
	
	VisuEngine->drawWindows();
}




void disturbPose()
{
	VectorXf poseNoise(6);
	poseNoise.setZero();
	srand (time(NULL));
	for(int j=0;j<nb_kf;j++)
	{
		for(int i=0;i<3;i++)poseNoise[i]+=0.02*((double)rand()/(double)RAND_MAX-0.5);
		for(int i=0;i<3;i++)poseNoise[i+3]+=0.1*((double)rand()/(double)RAND_MAX-0.5);		
		
		PoseEstim[j]=HomogeneousMatrix(poseNoise)* PoseEstim[j];
	}

}
void disturbDepth()
{
	srand (time(NULL));
	
	float min_depth=invDepthGT[0].invDepth;
	float max_depth=invDepthGT[0].invDepth;
	for(int i=0;i<invDepthGT.size();i++) 
	{
		if(invDepthGT[i].invDepth<min_depth)min_depth=invDepthGT[i].invDepth;
		if(invDepthGT[i].invDepth>max_depth)max_depth=invDepthGT[i].invDepth;
	}
	
	for(int i=0;i<invDepthGT.size();i++) 
	{
		invDepthEstim[i].invDepth+=0.1*(max_depth-min_depth)*((double)rand()/(double)RAND_MAX-0.5);
	}


}

int kf_disp=0;

void createRandomEstim()
{
	invDepthEstim.clear();
	
	//get depth min and max (optional just because up to scale)
	float min_depth=invDepthGT[0].invDepth;
	float max_depth=invDepthGT[0].invDepth;
	for(int i=0;i<invDepthGT.size();i++) 
	{
		if(invDepthGT[i].invDepth<min_depth)min_depth=invDepthGT[i].invDepth;
		if(invDepthGT[i].invDepth>max_depth)max_depth=invDepthGT[i].invDepth;
	}
	  
	srand (time(NULL));
	for(int i=0;i<invDepthGT.size();i++) 
	{
		PointInvDepth newEstim;
		newEstim.srcKf=0;
		newEstim.meterCoord=invDepthGT[i].meterCoord;
		newEstim.invDepth=min_depth+(max_depth-min_depth)*((double)rand()/(double)RAND_MAX);
		newEstim.invDepthCovar=0;
		invDepthEstim.push_back(newEstim);
	}
}

void EstimEqualGT()
{
	for(int i=0;i<nb_kf;i++)
		PoseEstim[i]=PoseGT[i];
	
	invDepthEstim.clear();
	for(int i=0;i<invDepthGT.size();i++) 
		invDepthEstim.push_back(invDepthGT[i]);
}

void keyboard_process(unsigned char key, int x, int y)
{
	switch(key) {
		case 'k':
		        kf_disp= (kf_disp+1)% nb_kf;
			char fileName[200];
			sprintf(fileName,"View from kf %d, ",kf_disp);
			InKfWindow->setTitle(fileName);
			break;

		case 'c'://esc
			createRandomEstim();
			break;
		case 'd'://esc
			disturbPose();
			break;
		case 'f'://esc
			disturbDepth();
			break;
		case 'e'://esc
			EstimEqualGT();
			break;
		case 'b'://esc
			DoMiniBA();
			break;			
		case 27://esc
			delete VisuEngine;
			exit(0);
			break;	
	  
	}

}

std::vector<Vector2f> noisyMeasures[nb_kf];

Vector2f ProjInvDepthPoint(PointInvDepth &pt,HomogeneousMatrix &pose)
{
	Vector3f coord_homog=toHomogeneous(pt.meterCoord);
	Vector3f rotCoordPlusTrans=pose.get_rotation()*coord_homog+pt.invDepth*pose.get_translation();
	return rotCoordPlusTrans.segment(0,2)/rotCoordPlusTrans[2];
}

void createMap()
{
// 	//create loop around point at distance z=dist from cam
	float angle=10;float dist=0.6;
	
	//create Points that will do measures
	int nb_point_all=2000;	
	
	for(int i=0;i<nb_point_all;i++)
	{
		//get coord
		Vector3f pt;
		float angle2=2*3.1415*((double)rand()/(double)RAND_MAX);
		pt[0]=cos(angle2)*(dist*2);
		pt[1]=0;
		pt[2]=sin(angle2)*(dist*2);
		
		for(int j=0;j<3;j++)pt[j]+=0.3*((double)rand()/(double)RAND_MAX-0.5);
		//std::cout<<((double)rand()/(double)RAND_MAX-0.5)<<std::endl;
		RandomPoints.push_back(pt);
	}	
	
	
	PoseRef.TranslateZ(dist);
	cv::Mat img(im_h,im_w,CV_8UC1);
	img.setTo(0);
	
	//create all Keyframes
	srand (time(NULL));
	for(int j=0;j<nb_kf;j++)
	{
		VectorXf poseNoise(6);poseNoise.setZero();
		for(int i=0;i<3;i++)poseNoise[i]+=0.1*((double)rand()/(double)RAND_MAX-0.5);
		for(int i=0;i<3;i++)poseNoise[i+3]+=0.15*((double)rand()/(double)RAND_MAX-0.5);		
		
		PoseGT[j]=HomogeneousMatrix(poseNoise)* PoseRef;
	}
	for(int j=0;j<nb_kf;j++)
		PoseEstim[j]=PoseRef;
	
	//create invDepth GT for each point that projects in ref
	for(int i=0;i<nb_point_all;i++)
	{
		Vector3f ptCam=PoseRef*RandomPoints[i];
		if(ptCam[2]>0)
		{
			Vector2f proj=myCamera.Project(ptCam);
			if(proj[0]>0 && proj[0]<im_w && proj[1]>0 && proj[1]<im_h)
			{
				PointInvDepth newInvPoint;
				newInvPoint.srcKf=0;
				newInvPoint.meterCoord=myCamera.ProjectZ1(ptCam);
				newInvPoint.invDepth=1./ptCam[2];
				newInvPoint.invDepthCovar=0;
				invDepthGT.push_back(newInvPoint);
			}
		}
	}
	
	//create noisy measures
	for(int j=0;j<nb_kf;j++)
	{
		HomogeneousMatrix relPoseGT=PoseGT[j]*PoseRef.inverse();
		for(int i=0;i<invDepthGT.size();i++)
		{
			Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT);//current projection
			noisyMeasures[j].push_back(x_d+getImagePlanNoise());
		}
	  
	}
	  
	  
	
	
}

//draw map viewed from kf_disp if disp_best_local ==0
//else draw map from best ministereo pair of kf_disp
void addDrawFunction(void) 
{	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glEnable(GL_POINT_SMOOTH);
	
	set2DGLProjection();
	glColor3f(0,1,0);
	

	glColor3f(1,1,1);
	unset2DGLProjection();
	//glutSwapBuffers();
}
struct miniBAjacz
{
      short index_pt_opt;//if max angle not big enough => pt not optim then index_pt_opt=-1
      short index_cam_opt;//if max angle not big enough => pt not optim then index_pt_opt=-1
      float weight;
      Vector2f proj_error;//deriv wrt point pos (only if index_pt_opt!=-1)
      Vector2f de_dz;//deriv wrt point pos (only if index_pt_opt!=-1)
      MatrixXf de_dp;//deriv wrt cam pos
};

void DoMiniBA()
{
	int nb_feat_ref=invDepthGT.size();
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;	
	
	HomogeneousMatrix relPoseGT[nb_kf];
	HomogeneousMatrix relPoseEstim[nb_kf];
	for(int k=0;k<nb_kf;k++)
	{
		relPoseGT[k]=PoseGT[k]*PoseRef.inverse();
		relPoseEstim[k]=PoseEstim[k]*PoseRef.inverse();
	}
	
	for(int iter=0;iter<1;iter++)
	{
		//get Tukey factor
		std::cout<<"iter "<<iter<<std::endl;
		//std::cout<<"get Tukey factor iter "<<iter<<std::endl;
		float sigma_tukey[nb_kf];
		for(int k=0;k<nb_kf;k++)
		{		
			std::vector<float> vdErrorSquared;
			for(int i=0;i<nb_feat_ref;i++)
			{

				//Vector3f mapPointGT=PoseGT[k]*(toHomogeneous(invDepthGT[i].meterCoord)/invDepthGT[i].invDepth);
				//Vector3f mapPointEstim=PoseEstim[k]*(toHomogeneous(invDepthEstim[i].meterCoord)/invDepthEstim[i].invDepth);
				//if(mapPointGT[2]>0 && mapPointEstim[2]>0)
				{
				    //Vector2f x_c=myCamera.ProjectZ1(mapPointEstim);//current projection
				    //Vector2f x_d=myCamera.ProjectZ1(mapPointGT);//current projection
				    //Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT[k]);//current projection
				    Vector2f x_d=noisyMeasures[k][i];
				    Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
				    Vector2f error=myCamera.m2PixProjJac()*(x_d-x_c);//error in pixels
				    vdErrorSquared.push_back(error.transpose()*error);
				}
				
			}
			//std::cout<<"nb meas view ["<<k<<"] = "<<vdErrorSquared.size()<<std::endl;
			sigma_tukey[k]=getSigmaSquared(vdErrorSquared);
			if(sigma_tukey[k] < MINSIGMATUKEY)
				sigma_tukey[k] = MINSIGMATUKEY;	
		}
		
		
		//now collect all jacobians for camera and depth features to optimise (for that needs to be matched, have a positive depth which is not outlier)
		//get LUT for depths to optimise
		//std::cout<<"get LUTs"<<std::endl;
		int LUT_to_opt[nb_feat_ref];
		int LUT_to_main[nb_feat_ref];
		
		bool first_fixedPoint_passed=false;
		int opt_cpt=0;
		for(int i=0;i<nb_feat_ref;i++)
		{
			LUT_to_opt[i]=-1;
			for(int k=0;k<nb_kf;k++)
			{
				//Vector3f mapPointGT=PoseGT[k]*(toHomogeneous(invDepthGT[i].meterCoord)/invDepthGT[i].invDepth);
				//Vector3f mapPointEstim=PoseEstim[k]*(toHomogeneous(invDepthEstim[i].meterCoord)/invDepthEstim[i].invDepth);
				//if(mapPointGT[2]>0 && mapPointEstim[2]>0)
				{
					//Vector2f x_c=myCamera.ProjectZ1(mapPointEstim);//current projection
					//Vector2f x_d=myCamera.ProjectZ1(mapPointGT);//current projection
					//Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT[k]);//current projection
					Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection

					Vector2f error=myCamera.m2PixProjJac()*(x_d-x_c);//error in pixels

					float TukeyCoef=squareRootTukey(error.transpose()*error,sigma_tukey[k]);
					if(TukeyCoef>0)
					{
						//check if variation of depth changes anything at all
						Vector3f coord_homog=toHomogeneous(invDepthEstim[i].meterCoord);
						Vector3f rotCoord=relPoseEstim[k].get_rotation()*coord_homog;
						//Vector3f rotCoordPlusTrans=relPoseEstim[k].get_rotation()*coord_homog+invDepthEstim[i].invDepth*relPoseEstim[k].get_translation();
						
						Vector3f t=relPoseEstim[k].get_translation();
						float w=invDepthEstim[i].invDepth;
						float Nx=rotCoord[0]+w*t[0];
						float Ny=rotCoord[1]+w*t[1];
						float D=1.+w*t[2];
						
						Vector2f de_dz;
						de_dz[0]=(-t[0]*D+t[2]*Nx)/(D*D);
						de_dz[1]=(-t[1]*D+t[2]*Ny)/(D*D);
						
						//if(de_dz[0]!=0 || de_dz[1]!=0)
						if(de_dz[0]*de_dz[0]>1e-10 || de_dz[1]*de_dz[1]>1e-10)
						{
							if(!first_fixedPoint_passed)//here choose the first one, night be better to choose one with largest de_dz...
								first_fixedPoint_passed=true;
							else
							{
								LUT_to_opt[i]=opt_cpt;
								LUT_to_main[opt_cpt]=i;
									
								opt_cpt++;
								break;//check next feature
							}
						}
					}

				}
			}
		}
		
		//get error and jacobians
		//std::cout<<"get error and jacs"<<std::endl;
		float residue_robust=0;
		float valid_wmeas=0;
		std::vector<miniBAjacz> fJacobian;		

		for(int k=0;k<nb_kf;k++)
		{		
			for(int i=0;i<nb_feat_ref;i++)
			{
				//Vector3f mapPointGT=relPoseGT[k]*(toHomogeneous(invDepthGT[i].meterCoord)/invDepthGT[i].invDepth);
				Vector3f mapPointEstim=relPoseEstim[k]*(toHomogeneous(invDepthEstim[i].meterCoord)/invDepthEstim[i].invDepth);
				//if(mapPointGT[2]>0 && mapPointEstim[2]>0)
				{
					//Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT[k]);//current projection
					Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
					Vector2f error=(x_d-x_c);//error in pixels
					Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels

					float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
					if(TukeyCoef>0)
					{
						residue_robust+=TukeyCoef*sqrt(errorPix.squaredNorm());
						valid_wmeas+=TukeyCoef;
							
						miniBAjacz newJc;
						newJc.proj_error=error;
						//newJc.weight=TukeyCoef;
						newJc.weight=1.;
						
						Vector3f coord_homog=toHomogeneous(invDepthEstim[i].meterCoord);
						Vector3f rotCoord=relPoseEstim[k].get_rotation()*coord_homog;
						//Vector3f rotCoordPlusTrans=relPoseEstim[k].get_rotation()*coord_homog+invDepthEstim[i].invDepth*relPoseEstim[k].get_translation();
						
						Vector3f t=relPoseEstim[k].get_translation();
						float w=invDepthEstim[i].invDepth;
						float Nx=rotCoord[0]+w*t[0];
						float Ny=rotCoord[1]+w*t[1];
						float D=1.+w*t[2];
						
						//compute jacobian with respect to cam position
						newJc.index_cam_opt=k;
						
						MatrixXf de_dp(2,6);
						de_dp(0,0)=-w/D;	de_dp(0,1)=0;		de_dp(0,2)=Nx*w/(D*D);	
						de_dp(1,0)=0;		de_dp(1,1)=-w/D;	de_dp(1,2)=Ny*w/(D*D);	
						
						de_dp(0,3)=rotCoord[1]*Nx/(D*D);	de_dp(0,4)=(-D-rotCoord[0]*Nx)/(D*D);	de_dp(0,5)=rotCoord[1]/D;
						de_dp(1,3)=(D+rotCoord[1]*Ny)/(D*D);	de_dp(1,4)=-rotCoord[0]*Ny/(D*D);	de_dp(1,5)=-rotCoord[1]/D;
						
						newJc.de_dp=de_dp;


						newJc.index_pt_opt=LUT_to_opt[i];
						//std::cout<<"newJc.index_pt_opt= "<<newJc.index_pt_opt <<std::endl;
						if(newJc.index_pt_opt!=-1)//hapen only with first point on which we fix depth so that we don t have up to scale problema
						{
							//get jacobien  of error with respect to vairation of depth
							Vector2f de_dz;
							de_dz[0]=(-t[0]*D+t[2]*Nx)/(D*D);
							de_dz[1]=(-t[1]*D+t[2]*Ny)/(D*D);
							newJc.de_dz=de_dz;

						}
						
						/*newJc.de_dp=-myCamera.ProjectZ1_Jac_Dp(mapPointEstim);

						newJc.index_pt_opt=LUT_to_opt[i];
						if(newJc.index_pt_opt!=-1)//hapen only with first point on which we fix depth so that we don t have up to scale problema
						{
							//get jacobien  of error with respect to vairation of  depth
							Vector3f JacXInvd;
							float inv_depth_cur=invDepthEstim[i].invDepth;
							JacXInvd[0]=-invDepthGT[i].meterCoord[0]/(inv_depth_cur*inv_depth_cur);
							JacXInvd[1]=-invDepthGT[i].meterCoord[1]/(inv_depth_cur*inv_depth_cur);
							JacXInvd[2]=-1./(inv_depth_cur*inv_depth_cur);
							newJc.de_dz=-myCamera.ProjectZ1_Jac_X(mapPointEstim)*(relPoseEstim[k]).get_rotation()*JacXInvd;	
						}*/
						
						fJacobian.push_back(newJc);

					}
				}
			}
	
		}
		std::cout<<"residue_robust = "<<residue_robust<<std::endl;
		std::cout<<"valid_wmeas = "<<valid_wmeas<<std::endl;
		
		//compute Hessain and update
		int nbCamsToUpdate=nb_kf;
		int nbPointsToUpdate=opt_cpt;
		VectorXf Jtex(nbPointsToUpdate);Jtex.setZero();
		VectorXf Jtep(6*nbCamsToUpdate);Jtep.setZero();

		float Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i]=0;
		MatrixXf Hpp(6*nbCamsToUpdate,6*nbCamsToUpdate);Hpp.setZero();
		MatrixXf Hxp(nbPointsToUpdate,6*nbCamsToUpdate);Hxp.setZero();
		
		
		for(int i=0;i<fJacobian.size();i++)			
		{
			short &pt_opt_id=fJacobian[i].index_pt_opt;
			short &pt_cam_id=fJacobian[i].index_cam_opt;

			if(pt_opt_id!=-1)
			{
				//update Jte
				Jtex[pt_opt_id]+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].proj_error;
				//update Hessian
				Hxx[pt_opt_id]+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].de_dz;
				
				Hxp.block(pt_opt_id,6*pt_cam_id,1,6)+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].de_dp;				
			}
			
			//update Jte
			Jtep.segment(6*pt_cam_id,6)+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].proj_error;
			//update Hessian
			Hpp.block(6*pt_cam_id,6*pt_cam_id,6,6)+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].de_dp;	
		}
		
		
		for(int i=0;i<nbPointsToUpdate;i++)
			Hxx[i]=(1.+mLMLambda)*Hxx[i];	
		for(int i=0;i<6*nbCamsToUpdate;i++)
			Hpp(i,i)=(1.+mLMLambda)*Hpp(i,i);	
		
		
		/*
		//without optim from Engel
		MatrixXf H(nbPointsToUpdate+6*nbCamsToUpdate,nbPointsToUpdate+6*nbCamsToUpdate);H.setZero();
		VectorXf J(nbPointsToUpdate+6*nbCamsToUpdate);
		
		for(int i=0;i<nbPointsToUpdate;i++)J[i]=Jtex[i];
		for(int i=0;i<6*nbCamsToUpdate;i++)J[nbPointsToUpdate+i]=Jtep[i];
	
		for(int k=0;k<nbPointsToUpdate;k++)
			H(k,k)=Hxx[k];
				
		for(int i=0;i<6*nbCamsToUpdate;i++)
			for(int j=0;j<6*nbCamsToUpdate;j++)H(nbPointsToUpdate+i,nbPointsToUpdate+j)=Hpp(i,j);
			
		for(int i=0;i<nbPointsToUpdate;i++)
			for(int j=0;j<6*nbCamsToUpdate;j++)
			{
				H(i,nbPointsToUpdate+j)=Hxp(i,j);
				H(nbPointsToUpdate+j,i)=Hxp(i,j);
			}
		
		//compute update
		Eigen::FullPivLU<MatrixXf> lu(H);
		VectorXf Dv(nbPointsToUpdate+6*nbCamsToUpdate);
		Dv=-(lu.inverse()*J);

		//get update cam positions
		VectorXf Dc(6*nbCamsToUpdate);
		Dc=Dv.segment(nbPointsToUpdate,6*nbCamsToUpdate);
		
		
		//get update inv depth positions
		VectorXf Dz_all(nbPointsToUpdate);
		Dz_all=Dv.segment(0,nbPointsToUpdate);
		
		*/

		//with optim from Engel
		float HxxInv[nbPointsToUpdate];
		for(int i=0;i<nbPointsToUpdate;i++)HxxInv[i]=1./Hxx[i];
		
		MatrixXf HxxInvHxp(nbPointsToUpdate,6*nbCamsToUpdate);//HxxInvHxp=HxxInv*Hxp;
		for(int i=0;i<nbPointsToUpdate;i++)
			for(int c=0;c<nbCamsToUpdate;c++)HxxInvHxp.block(i,6*c,1,6)=HxxInv[i]*Hxp.block(i,6*c,1,6);
		
		VectorXf HxxInvJtp(nbPointsToUpdate);//HxxInvJtp=HxxInv*Jtex;
		for(int i=0;i<nbPointsToUpdate;i++)HxxInvJtp.segment(i,1)=HxxInv[i]*Jtex.segment(i,1);
		
		MatrixXf A(6*nbCamsToUpdate,6*nbCamsToUpdate);	A=Hpp-Hxp.transpose()*HxxInvHxp;
		VectorXf B(6*nbCamsToUpdate);		B=Jtep-Hxp.transpose()*HxxInvJtp;

		Eigen::FullPivLU<MatrixXf> luAc(A);
		VectorXf Dc=-luAc.inverse()*B;	
		
		//update point coord and depths
		VectorXf Dz_all(nbPointsToUpdate);
		Dz_all=-HxxInvJtp-HxxInvHxp*Dc;
		
		HomogeneousMatrix newRelPoseEstim[nb_kf];
		for(int k=0;k<nb_kf;k++)
		{
			VectorXf Dci=0.2*Dc.segment(6*k,6);
			//std::cout<<"Dci["<<k<<"] = "<<Dci<<std::endl;
			//newRelPoseEstim[k]=HomogeneousMatrix(Dci)*relPoseEstim[k];
			//use rotation in compositional and translkation in additional
			newRelPoseEstim[k].set_rotation(HomogeneousMatrix(Dci).get_rotation()*relPoseEstim[k].get_rotation());
			newRelPoseEstim[k].set_translation(HomogeneousMatrix(Dci).get_translation()+relPoseEstim[k].get_translation());
		}

		//update point coord and depths
		PointInvDepth newInvDepthEstim[nb_feat_ref];
		for(int k=0;k<nb_feat_ref;k++)newInvDepthEstim[k]=invDepthEstim[k];

		for(int id_opt=0;id_opt<nbPointsToUpdate;id_opt++)
		{
			short id_glob=LUT_to_main[id_opt];
			float Dz=0.2*Dz_all[id_opt];
			//std::cout<<"Dz["<<id_glob<<"] = "<<Dz<<std::endl;
			newInvDepthEstim[id_glob].invDepth+=Dz;
			if(newInvDepthEstim[id_glob].invDepth<0)newInvDepthEstim[id_glob].invDepth=invDepthEstim[id_glob].invDepth;
			if(newInvDepthEstim[id_glob].invDepth<1e-6)newInvDepthEstim[id_glob].invDepth=1e-6;
			
		}
			
		
		
		/*float residue_robustAfter=0;
		float valid_wmeasAfter=0;
		for(int k=0;k<nb_kf;k++)
		{		
			for(int i=0;i<nb_feat_ref;i++)
			{
				//Vector3f mapPointGT=PoseGT[k]*(toHomogeneous(invDepthGT[i].meterCoord)/invDepthGT[i].invDepth);
				//Vector3f mapPointEstim=newPoseEstim[k]*(toHomogeneous(invDepthEstim[i].meterCoord)/invDepthEstim[i].invDepth);
				//if(mapPointGT[2]>0 && mapPointEstim[2]>0)
				{
					//Vector2f x_c=myCamera.ProjectZ1(mapPointEstim);//current projection
					//Vector2f x_d=myCamera.ProjectZ1(mapPointGT);//current projection
					//Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT[k]);//current projection
					Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_c=ProjInvDepthPoint(newInvDepthEstim[i],newRelPoseEstim[k]);//current projection
					Vector2f error=(x_d-x_c);//error in pixels
					Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels

					float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
					if(TukeyCoef>0)
					{
						residue_robustAfter+=TukeyCoef*sqrt(errorPix.squaredNorm());
						valid_wmeasAfter+=TukeyCoef;
					}

				}
			}
		}
		std::cout<<"residue_robustAfter = "<<residue_robust<<std::endl;
		std::cout<<"valid_wmeasAfter = "<<valid_wmeas<<std::endl;
		
		if(valid_wmeasAfter>0 && residue_robustAfter/valid_wmeasAfter<residue_robust/valid_wmeas)
		{*/
			for(int k=0;k<nb_feat_ref;k++)invDepthEstim[k]=newInvDepthEstim[k];
			for(int k=0;k<nb_kf;k++)relPoseEstim[k]=newRelPoseEstim[k];
			
		/*	mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
		}
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}*/
	}
	//rescale with mean invDepth
	float meanInvDepthGT=0;
	for(int k=0;k<nb_feat_ref;k++)meanInvDepthGT=invDepthGT[k].invDepth;
	float meanInvDepthEstim=0;
	for(int k=0;k<nb_feat_ref;k++)meanInvDepthEstim=invDepthEstim[k].invDepth;
	
	for(int k=0;k<nb_feat_ref;k++)invDepthEstim[k].invDepth=invDepthEstim[k].invDepth*meanInvDepthGT/meanInvDepthEstim;
	for(int k=0;k<nb_kf;k++)relPoseEstim[k].set_translation(relPoseEstim[k].get_translation()*meanInvDepthEstim/meanInvDepthGT);
	
	
	for(int k=0;k<nb_kf;k++)
		PoseEstim[k]=relPoseEstim[k]*PoseRef;
	
	
	
	
}
