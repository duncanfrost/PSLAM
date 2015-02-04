//perform a live BA directly from small movement of camera and tracking of KLT features
//see [3D Reconstruction from Accidental Motion - Fisher Yu]

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"

#include <opencv2/video/tracking.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;

void Idle(void) ;
void processNormalKeys(unsigned char key, int x, int y);
void addDrawFunction(void) ;
void useNewMatches(std::vector<p_match> &_matchesCurrent);
void EstimateCurrentPose(std::vector<p_match> &_matchesCurrent,HomogeneousMatrix &relPoseBA);
std::vector<p_match> matchFeatures(cv::Mat &_img);

//Visualization
VisualisationModule *VisuEngine;
MapWindow *mMapWindow;
//image acquisition
VideoSourceLiveCV *myVideoSource;
Camera myCamera;//camera object, calibration ...

//reference
HomogeneousMatrix PoseRef;
	

//reference
cv::Mat imgRef;	
cv::Mat imgPrev;

//tracking info
vector<Point2f> points[2];
Size subPixWinSize(10,10), winSize(31,31);
TermCriteria termcrit(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 20, 0.03);
const int MAX_COUNT = 500;
bool nightMode = false;
std::vector<p_match> matchesCurrent;

void initMatcher();

//all matches per feature
Vector2f **AllMatches;

//reslut from BA
#define nb_kf_max 20
int nb_feat_ref;
float min_idepth,max_idepth;
std::vector<PointInvDepth> invDepthEstim;
int nb_kf=0;
HomogeneousMatrix PoseEstim[nb_kf_max];
int DoMiniBA(std::vector<int> &fix_kf,std::vector<int> &fix_pt);
int optimCamPose(std::vector<int> &ptToUse);
int optimPointPose();

HomogeneousMatrix poseEstimFromMap;

Vector2f ProjInvDepthPoint(PointInvDepth &pt,HomogeneousMatrix &pose)
{
	Vector3f coord_homog=toHomogeneous(pt.meterCoord);
	Vector3f rotCoordPlusTrans=pose.get_rotation()*coord_homog+pt.invDepth*pose.get_translation();
	return rotCoordPlusTrans.segment(0,2)/rotCoordPlusTrans[2];
}

int main(int argc, char** argv)
{

    #ifdef USE_OMP_C
        coutGreen << "We are using openMP" << endlGreen;
    #endif


	int im_w=640;
	int im_h=480;
	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource->initCam();
	myCamera=*myVideoSource->getPointerCamera();

	//get reference image
	myVideoSource->grabNewFrame();
	myVideoSource->GetFramePointerBW()->copyTo(imgRef);
	initMatcher();
	
	//init inv depth features
	srand (time(NULL));
	min_idepth=0.01;
	max_idepth=1.;
	nb_feat_ref=points[0].size();
	for(int i=0;i<nb_feat_ref;i++)
	{
		PointInvDepth newInvDepth;
		newInvDepth.srcKf=0;
		newInvDepth.meterCoord=myCamera.ToMeters(Vector2f(points[0][i].x,points[0][i].y));
		newInvDepth.invDepth=min_idepth+(max_idepth-min_idepth)*((double)rand()/(double)RAND_MAX);
		invDepthEstim.push_back(newInvDepth);
	}
	
	//alloc AllMatches
	AllMatches= new Vector2f*[nb_feat_ref];
	for(int i=0;i<nb_feat_ref;i++)AllMatches[i] = new Vector2f[nb_kf_max];
	
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeys);


	VisuEngine->addWindowMap("Map",640,480,&myCamera,mMapWindow);
	VisuEngine->setOnKeyPress(processNormalKeys);
	
	//display gt map and estimated map
	mMapWindow->addCamera(&PoseRef,Vector3f(1,0,0),2);
	for(int i=0;i<nb_kf_max;i++)
	{
		mMapWindow->addCamera(&PoseEstim[i],Vector3f(0,1,0),1);
	}
	mMapWindow->addCamera(&poseEstimFromMap,Vector3f(0,0,1),2);
	
	mMapWindow->addPointCloud(&invDepthEstim);
	
	//HomogeneousMatrix moveCam(0,0.8,0.8,0.7,0,0);
	//HomogeneousMatrix moveCam(0,2.1,2.1,2.0,0,0);
	HomogeneousMatrix moveCam;
	moveCam.TranslateZ(-0.3);
	moveCam.TranslateY(-0.5);
	moveCam.RotateX(-45);
	mMapWindow->moveCamera(moveCam);

	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}

bool poseProcess=false;
void Idle(void) 
{
	myVideoSource->grabNewFrame();
	cv::Mat &current_img_BW=*myVideoSource->GetFramePointerBW();
	
	if(!poseProcess && nb_kf<nb_kf_max)
	{
		std::cout<<"New image"<<std::endl;
		matchesCurrent=matchFeatures(current_img_BW);
		//also set matchesCurrentBucket
		//std::cout<<"Homography = "<<Homography<<std::endl;
		std::cout<<"matchesCurrent.size() = "<<matchesCurrent.size()<<std::endl;
		
		useNewMatches(matchesCurrent);
	}
	
	VisuEngine->drawWindows();
}

void useNewMatches(std::vector<p_match> &_matchesCurrent)
{
	//update AllMatches
	if(nb_kf<nb_kf_max)
	{
		for(int i=0;i<nb_feat_ref;i++)AllMatches[i][nb_kf]=Vector2f(-1,-1);
		for(int i=0;i<_matchesCurrent.size();i++)
			AllMatches[_matchesCurrent[i].i1p][nb_kf]=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
		nb_kf++;
	
		//check if algebraic decomposition of fundamental gives better reprojection error
		HomogeneousMatrix relPoseAlgeb;
		//bool isMotionGood=CheckNewStereo(_matchesCurrentBucket,&myCamera,relPoseAlgeb);
		bool isMotionGood=CheckNewStereo(_matchesCurrent,&myCamera,relPoseAlgeb);
		if(isMotionGood)
		{
			std::cout<<"motion is good"<<std::endl;
			relPoseAlgeb.set_translation(0.1*relPoseAlgeb.get_translation());
			
			//get pose current from min reprojection
			HomogeneousMatrix relPoseBA;
			EstimateCurrentPose(_matchesCurrent,relPoseBA);
			poseEstimFromMap=relPoseBA;

		  
			//compute reprojection error
			float reproj_error_BA=0;
			float reproj_error_Alg=0;
			
			float min_idepth_Alg=-1;
			float max_idepth_Alg=-1;
			for(int i=0;i<_matchesCurrent.size();i++)
			{
				Vector2f mes_p=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1p,_matchesCurrent[i].v1p));
				Vector2f mes_c=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
				float invdepth;
				float recAngle;
				int rec=reconstructionFromRaysInvDepth(mes_p,mes_c,relPoseAlgeb,invdepth, recAngle);//return 0 if points at infinity
				
				if(rec!=-1)
				{
					Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
					PointInvDepth currentInvPointAlg;
					currentInvPointAlg=invDepthEstim[_matchesCurrent[i].i1p];
					currentInvPointAlg.invDepth=invdepth;
					
					Vector2f x_c_alg=ProjInvDepthPoint(currentInvPointAlg,relPoseAlgeb);//current projection
					Vector2f x_c_BA=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[i].i1p],relPoseBA);//current projection
					
					Vector2f error_alg=x_d-x_c_alg;
					Vector2f error_BA=x_d-x_c_BA;
					
					reproj_error_Alg+=error_alg.transpose()*error_alg;
					reproj_error_BA+=error_BA.transpose()*error_BA;
					
					if(min_idepth_Alg==-1 || min_idepth_Alg>invdepth)min_idepth_Alg=invdepth;
					if(max_idepth_Alg==-1 || max_idepth_Alg<invdepth)max_idepth_Alg=invdepth;
				}
			}
			std::cout<<"reproj_error_Alg = "<<reproj_error_Alg<<std::endl;
			std::cout<<"reproj_error_BA = "<<reproj_error_BA<<std::endl;
			
			if(reproj_error_Alg<reproj_error_BA)
			{
				coutRed<<"Algebraic solution better"<<endlRed;
				//reset poses to identity
				for(int i=0;i<nb_kf;i++)PoseEstim[i].SetIdentity();
				PoseEstim[nb_kf-1]=relPoseAlgeb;
				
				//reset all features to random depth
				for(int i=0;i<nb_feat_ref;i++)
					invDepthEstim[i].invDepth=0;
					//invDepthEstim[i].invDepth=min_idepth_Alg+(max_idepth_Alg-min_idepth_Alg)*((double)rand()/(double)RAND_MAX);
				
				for(int i=0;i<_matchesCurrent.size();i++)
				{
					Vector2f mes_p=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1p,_matchesCurrent[i].v1p));
					Vector2f mes_c=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
					float invdepth;
					float recAngle;
					int rec=reconstructionFromRaysInvDepth(mes_p,mes_c,relPoseAlgeb,invdepth, recAngle);//return 0 if points at infinity
					
					if(rec!=-1)
						invDepthEstim[_matchesCurrent[i].i1p].invDepth=invdepth;
				}
				
				float reproj_error_test=0;
				for(int i=0;i<_matchesCurrent.size();i++)
				{
						Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
						Vector2f x_c_BA=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[i].i1p],PoseEstim[nb_kf-1]);//current projection
						if(i==0)std::cout<<"invDepthEstim[_matchesCurrent[i].i1p].invDepth = "<<invDepthEstim[_matchesCurrent[i].i1p].invDepth<<std::endl;
						
						Vector2f error_BA=x_d-x_c_BA;
						reproj_error_test+=error_BA.transpose()*error_BA;
				}
				std::cout<<"reproj_error_test = "<<reproj_error_test<<std::endl;
				
				
				//do BA fixing current KF and points
				std::vector<int> fix_kf;
				fix_kf.push_back(nb_kf-1);
				std::vector<int> fix_pt;
				for(int i=0;i<_matchesCurrent.size();i++)fix_pt.push_back(_matchesCurrent[i].i1p);
				//DoMiniBA(fix_kf,fix_pt);
				optimCamPose(fix_pt);
				//optimPointPose();

				reproj_error_test=0;
				for(int i=0;i<_matchesCurrent.size();i++)
				{

						Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));
						Vector2f x_c_BA=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[i].i1p],PoseEstim[nb_kf-1]);//current projection
						if(i==0)std::cout<<"invDepthEstim[_matchesCurrent[i].i1p].invDepth = "<<invDepthEstim[_matchesCurrent[i].i1p].invDepth<<std::endl;
						
						Vector2f error_BA=x_d-x_c_BA;
						reproj_error_test+=error_BA.transpose()*error_BA;
				}
				//shoud be the same since points and cam fixed in BA
				std::cout<<"reproj_error_test = "<<reproj_error_test<<std::endl;
				
			  
			}
		}
	}
	//do BA
	std::vector<int> fix_kf;
	std::vector<int> fix_pt;
	//DoMiniBA(fix_kf,fix_pt);
	
}

void processNormalKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case 27://esc
			exit(0);
			for(int i=0;i<nb_feat_ref;i++)delete[] AllMatches[i];
			delete[] AllMatches;
			break;
		case ' '://esc
			poseProcess=!poseProcess;
			break;
	}


	glutPostRedisplay();
}

void initMatcher()
{
	goodFeaturesToTrack(imgRef, points[0], MAX_COUNT, 0.01, 10, Mat(), 3, 0, 0.04);
	cornerSubPix(imgRef, points[0], subPixWinSize, Size(-1,-1), termcrit);
	imgRef.copyTo(imgPrev);
	points[1]=points[0];
}

std::vector<p_match> matchFeatures(cv::Mat &_img)
{
	vector<uchar> status;
        vector<float> err;
	
	vector<Point2f> pointsTemp=points[1];
	calcOpticalFlowPyrLK(imgPrev, _img, pointsTemp, points[1], status, err, winSize,
                                 3, termcrit, 0, 0.001);
	_img.copyTo(imgPrev);

	 std::vector<p_match> matches;
	for(int i=0;i<points[1].size();i++)
	{
		p_match newMatch(points[0][i].x,points[0][i].y,i,points[1][i].x,points[1][i].y,i);
		matches.push_back(newMatch);
	}

	return matches;
}

void addDrawFunction(void) 
{	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	set2DGLProjection();
	glPointSize(3.0);
	glColor3f(0,1,0);

		
	for(int i=0;i<matchesCurrent.size();i++)
	{
		glLineWidth(1.0);
		glColor3f(0,0.,1);
		glBegin(GL_LINES);	
		glVertex2f(matchesCurrent[i].u1p,matchesCurrent[i].v1p);
		glColor3f(0,1.,1);
		glVertex2f(matchesCurrent[i].u1c,matchesCurrent[i].v1c);
		glEnd();
	}

	
	glColor3f(1,1,1);
	unset2DGLProjection();

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
int DoMiniBA(std::vector<int> &fix_kf,std::vector<int> &fix_pt)
{
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;	
	
	HomogeneousMatrix relPoseEstim[nb_kf];
	for(int k=0;k<nb_kf;k++)
	{
		relPoseEstim[k]=PoseEstim[k]*PoseRef.inverse();
	}
	
	for(int iter=0;iter<10;iter++)
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
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f error=myCamera.m2PixProjJac()*(x_d-x_c);//error in pixels
						vdErrorSquared.push_back(error.transpose()*error);
					}
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
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
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
									if(std::find(fix_pt.begin(),fix_pt.end(),i)==fix_pt.end())
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
			}
		}
		int LUTcam_to_opt[nb_kf];
		int LUTcam_to_main[nb_kf];
		int cam_opt_cpt=0;
		for(int k=0;k<nb_kf;k++)
		{
			if(std::find(fix_kf.begin(),fix_kf.end(),k)==fix_kf.end())
			{
				LUTcam_to_opt[k]=cam_opt_cpt;
				LUTcam_to_main[cam_opt_cpt]=k;
				cam_opt_cpt++;
			}
			else
				LUTcam_to_opt[k]=-1;
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
					//Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
						Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels

						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
							//residue_robust+=TukeyCoef*sqrt(errorPix.squaredNorm());
							//valid_wmeas+=TukeyCoef;
							residue_robust+=error.transpose()*error;
							valid_wmeas++;
								
							miniBAjacz newJc;
							newJc.proj_error=error;
							//newJc.weight=TukeyCoef;
							newJc.weight=1.;
							if(std::find(fix_pt.begin(),fix_pt.end(),i)==fix_pt.end())
								newJc.weight=10;
							
							Vector3f coord_homog=toHomogeneous(invDepthEstim[i].meterCoord);
							Vector3f rotCoord=relPoseEstim[k].get_rotation()*coord_homog;
							//Vector3f rotCoordPlusTrans=relPoseEstim[k].get_rotation()*coord_homog+invDepthEstim[i].invDepth*relPoseEstim[k].get_translation();
							
							Vector3f t=relPoseEstim[k].get_translation();
							float w=invDepthEstim[i].invDepth;
							float Nx=rotCoord[0]+w*t[0];
							float Ny=rotCoord[1]+w*t[1];
							float D=1.+w*t[2];
							
							//compute jacobian with respect to cam position
							newJc.index_cam_opt=LUTcam_to_opt[k];
							if(newJc.index_cam_opt!=-1)
							{
							
								MatrixXf de_dp(2,6);
								de_dp(0,0)=-w/D;	de_dp(0,1)=0;		de_dp(0,2)=Nx*w/(D*D);	
								de_dp(1,0)=0;		de_dp(1,1)=-w/D;	de_dp(1,2)=Ny*w/(D*D);	
								
								de_dp(0,3)=rotCoord[1]*Nx/(D*D);	de_dp(0,4)=(-D-rotCoord[0]*Nx)/(D*D);	de_dp(0,5)=rotCoord[1]/D;
								de_dp(1,3)=(D+rotCoord[1]*Ny)/(D*D);	de_dp(1,4)=-rotCoord[0]*Ny/(D*D);	de_dp(1,5)=-rotCoord[1]/D;
								
								newJc.de_dp=de_dp;
							}

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
						else
						{
							residue_robust+=sigma_tukey[k];
							valid_wmeas++;
						}
					}
				}
			}
	
		}
		std::cout<<"residue_robust = "<<residue_robust<<std::endl;
		std::cout<<"valid_wmeas = "<<valid_wmeas<<std::endl;
		
		//compute Hessain and update
		int nbCamsToUpdate=nb_kf-fix_kf.size();
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
				if(pt_cam_id!=-1)
					Hxp.block(pt_opt_id,6*pt_cam_id,1,6)+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].de_dp;				
			}
			if(pt_cam_id!=-1)
			{
				//update Jte
				Jtep.segment(6*pt_cam_id,6)+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].proj_error;
				//update Hessian
				Hpp.block(6*pt_cam_id,6*pt_cam_id,6,6)+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].de_dp;	
			}
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
		for(int k=0;k<nb_kf;k++)newRelPoseEstim[k]=relPoseEstim[k];
		  
		for(int k=0;k<nbCamsToUpdate;k++)
		{
			VectorXf Dci=0.2*Dc.segment(6*k,6);
			if(isnan(Dci[0])){coutRed<<"nanerie cam"<<endlRed;return -1;};
			short id_glob=LUTcam_to_main[k];
			//std::cout<<"Dci["<<k<<"] = "<<Dci<<std::endl;
			//newRelPoseEstim[k]=HomogeneousMatrix(Dci)*relPoseEstim[k];
			//use rotation in compositional and translkation in additional
			newRelPoseEstim[id_glob].set_rotation(HomogeneousMatrix(Dci).get_rotation()*relPoseEstim[id_glob].get_rotation());
			newRelPoseEstim[id_glob].set_translation(HomogeneousMatrix(Dci).get_translation()+relPoseEstim[id_glob].get_translation());
		}

		//update point coord and depths
		PointInvDepth newInvDepthEstim[nb_feat_ref];
		for(int k=0;k<nb_feat_ref;k++)newInvDepthEstim[k]=invDepthEstim[k];

		for(int id_opt=0;id_opt<nbPointsToUpdate;id_opt++)
		{
			short id_glob=LUT_to_main[id_opt];
			float Dz=0.2*Dz_all[id_opt];
			if(isnan(Dz) || isinf(Dz)){coutRed<<"nanerie iz"<<endlRed;return -1;};
			//std::cout<<"Dz["<<id_glob<<"] = "<<Dz<<std::endl;
			newInvDepthEstim[id_glob].invDepth+=Dz;
			if(newInvDepthEstim[id_glob].invDepth<0)newInvDepthEstim[id_glob].invDepth=invDepthEstim[id_glob].invDepth;
			if(newInvDepthEstim[id_glob].invDepth<1e-6)newInvDepthEstim[id_glob].invDepth=1e-6;
			if(isnan(newInvDepthEstim[id_glob].invDepth) || isinf(newInvDepthEstim[id_glob].invDepth)){coutRed<<"nanerie iz"<<endlRed;return -1;};
			
			
		}
			
		
		
		float residue_robustAfter=0;
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
					//Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(newInvDepthEstim[i],newRelPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
						Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels

						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
							//residue_robustAfter+=TukeyCoef*sqrt(errorPix.squaredNorm());
							//valid_wmeasAfter+=TukeyCoef;
							residue_robustAfter+=error.transpose()*error;
							valid_wmeasAfter++;
						}
						else
						{
							residue_robustAfter+=sigma_tukey[k];
							valid_wmeasAfter++;
						}
					}

				}
			}
		}
		std::cout<<"residue_robustAfter = "<<residue_robust<<std::endl;
		std::cout<<"valid_wmeasAfter = "<<valid_wmeas<<std::endl;
		
		if(valid_wmeasAfter>0 && residue_robustAfter/valid_wmeasAfter<residue_robust/valid_wmeas)
		{
			for(int k=0;k<nb_feat_ref;k++)invDepthEstim[k]=newInvDepthEstim[k];
			for(int k=0;k<nb_kf;k++)relPoseEstim[k]=newRelPoseEstim[k];
			
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
		}
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}
	}
	//rescale with mean invDepth
	float meanInvDepthGT=1.;
	float meanInvDepthEstim=0;
	for(int k=0;k<nb_feat_ref;k++)meanInvDepthEstim=invDepthEstim[k].invDepth;
	//std::cout<<"meanInvDepthEstim = "<<meanInvDepthEstim<<std::endl;
	
	if(meanInvDepthEstim!=0)
	{
	for(int k=0;k<nb_feat_ref;k++)invDepthEstim[k].invDepth=invDepthEstim[k].invDepth*meanInvDepthGT/meanInvDepthEstim;
	for(int k=0;k<nb_kf;k++)relPoseEstim[k].set_translation(relPoseEstim[k].get_translation()*meanInvDepthEstim/meanInvDepthGT);
	}
	
	
	for(int k=0;k<nb_kf;k++)
		PoseEstim[k]=relPoseEstim[k]*PoseRef;
}

void EstimateCurrentPose(std::vector<p_match> &_matchesCurrent,HomogeneousMatrix &relPose)
{
	int nb_iter=30;
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;
	
	for(int iter=0;iter<nb_iter;iter++)
	{
		//get tukey factor for robust estimation
		std::vector<float> vdErrorSquared;
		for(int m=0;m<_matchesCurrent.size();m++)
		{
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[m].u1c,_matchesCurrent[m].v1c));//desired projection= measurement
			Vector2f x_c=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[m].i1p],relPose);//current projection
			
			Vector2f errorPix=myCamera.m2PixProjJac()*(x_d-x_c);
			//if(measure.fromPoint==0)std::cout<<"\tLocal: "<<myCamera.ToPixels(x_d).transpose()<<"  \tcurrent proj: "<<myCamera.ToPixels(x_c).transpose()<<std::endl;
			//if(measure.fromPoint==1)std::cout<<"\tPoint: "<<myCamera.ToPixels(x_d).transpose()<<"  \tcurrent proj: "<<myCamera.ToPixels(x_c).transpose()<<std::endl;
			vdErrorSquared.push_back(errorPix.transpose()*errorPix);
			
		}

		if(vdErrorSquared.size()==0)
		{
			coutRed<<"no measure for pose computation"<<endlRed;
			break;
		}
		
		float sigma_tukey=getSigmaSquared(vdErrorSquared);
		if(sigma_tukey < MINSIGMATUKEY)
			sigma_tukey = MINSIGMATUKEY;
		
		//std::cout<<"\tsigma_tukey = "<<sqrt(sigma_tukey)<<std::endl;
	  
	  
		VectorXf Jte(6);Jte.setZero();
		MatrixXf H(6,6);H.setZero();
		
		float residue=0;
		float weightTotal=0;
		int nb_meas_used=0;
		for(int i=0;i<_matchesCurrent.size();i++)
		{
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));//desired projection= measurement
			Vector2f x_c=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[i].i1p],relPose);//current projection
			
			Vector2f error=x_d-x_c;
			float norm_reproj_error=(myCamera.m2PixProjJac()*error).squaredNorm();
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);					
			if(TukeyCoef>0 && invDepthEstim[_matchesCurrent[i].i1p].invDepth!=0)
			{
				float norm_reproj_error=error.transpose()*error;						

				//get jacobien  of error with respect to variation of camera pose
				Vector3f mapPointsCam=relPose*(toHomogeneous(invDepthEstim[_matchesCurrent[i].i1p].meterCoord)/invDepthEstim[_matchesCurrent[i].i1p].invDepth);
				//MatrixXf de_dp=-myCamera.ProjectZ1_Jac_Dp(mapPointsCam);
				
				Vector3f coord_homog=toHomogeneous(invDepthEstim[_matchesCurrent[i].i1p].meterCoord);
				Vector3f rotCoord=relPose.get_rotation()*coord_homog;
							//Vector3f rotCoordPlusTrans=relPoseEstim[k].get_rotation()*coord_homog+invDepthEstim[i].invDepth*relPoseEstim[k].get_translation();
							
				Vector3f t=relPose.get_translation();
				float w=invDepthEstim[_matchesCurrent[i].i1p].invDepth;
				float Nx=rotCoord[0]+w*t[0];
				float Ny=rotCoord[1]+w*t[1];
				float D=1.+w*t[2];

				MatrixXf de_dp(2,6);
				de_dp(0,0)=-w/D;	de_dp(0,1)=0;		de_dp(0,2)=Nx*w/(D*D);	
				de_dp(1,0)=0;		de_dp(1,1)=-w/D;	de_dp(1,2)=Ny*w/(D*D);	
				
				de_dp(0,3)=rotCoord[1]*Nx/(D*D);	de_dp(0,4)=(-D-rotCoord[0]*Nx)/(D*D);	de_dp(0,5)=rotCoord[1]/D;
				de_dp(1,3)=(D+rotCoord[1]*Ny)/(D*D);	de_dp(1,4)=-rotCoord[0]*Ny/(D*D);	de_dp(1,5)=-rotCoord[1]/D;

				
				
				Jte+=de_dp.transpose()*error;
				
				H+=de_dp.transpose()*de_dp;
				
				float err_pix=sqrt((myCamera.m2PixProjJac()*error).squaredNorm());
				residue+=TukeyCoef*err_pix;
				weightTotal+=TukeyCoef;
				nb_meas_used++;
			}
					
		}

		if(nb_meas_used<4)break;
		
		for(int i=0;i<6;i++)
			H(i,i)=(1.+mLMLambda)*H(i,i);
		
		Eigen::FullPivLU<MatrixXf> lu(H);
		//todo try p+=Dp
		VectorXf Dp(6);
		Dp=-0.3*(lu.inverse()*Jte);
		if(!isnan(Dp[0]))
		{
		
			HomogeneousMatrix relPoseNew=HomogeneousMatrix(Dp)* relPose;
			
			float residueAfter=0;
			float weightTotalAfter=0;
			for(int i=0;i<_matchesCurrent.size();i++)
			{
				//matrix to be filled			
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=myCamera.ToMeters(Vector2f(_matchesCurrent[i].u1c,_matchesCurrent[i].v1c));//desired projection= measurement
				Vector2f x_c=ProjInvDepthPoint(invDepthEstim[_matchesCurrent[i].i1p],relPoseNew);//current projection
				
				Vector2f error=x_d-x_c;
				float norm_reproj_error=(myCamera.m2PixProjJac()*error).squaredNorm();
				float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);					
				if(TukeyCoef>0 && invDepthEstim[_matchesCurrent[i].i1p].invDepth!=0)
				{
					float err_pix=sqrt((myCamera.m2PixProjJac()*error).squaredNorm());
					weightTotalAfter+=TukeyCoef;
					residueAfter+=TukeyCoef*err_pix;
				}
				
				
			}
			//std::cout<<"\tresidue/weightTotal "<<residue/weightTotal<<std::endl;
			//std::cout<<"\tresidueAfter/weightTotalAfter "<<residueAfter/weightTotalAfter<<std::endl;
			if(weightTotalAfter>0 && residueAfter/weightTotalAfter<residue/weightTotal)
			{
				mdLambdaFactor = 2.0;
				mLMLambda *= 0.3;
				relPose=relPoseNew;
			}
			else
			{
				mLMLambda = mLMLambda * mdLambdaFactor;
				mdLambdaFactor = mdLambdaFactor * 2;
			}
		
		}
		else break;
		
		
	}  
}
int optimCamPose(std::vector<int> &ptToUse)
{
	std::cout<<"optimCamPose"<<std::endl;
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;	
	
	HomogeneousMatrix relPoseEstim[nb_kf];
	for(int k=0;k<nb_kf;k++)
	{
		relPoseEstim[k]=PoseEstim[k]*PoseRef.inverse();
	}
	
	for(int iter=0;iter<10;iter++)
	{
		//get Tukey factor
		std::cout<<"iter "<<iter<<std::endl;
		//std::cout<<"get Tukey factor iter "<<iter<<std::endl;
		float sigma_tukey[nb_kf];
		int nb_meas_per_KF[nb_kf];
		for(int k=0;k<nb_kf;k++)
		{		
			nb_meas_per_KF[k]=0;
			std::vector<float> vdErrorSquared;
			for(int i=0;i<nb_feat_ref;i++)
				if(std::find(ptToUse.begin(),ptToUse.end(),i)!=ptToUse.end())
			{

				{
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f error=myCamera.m2PixProjJac()*(x_d-x_c);//error in pixels
						vdErrorSquared.push_back(error.transpose()*error);
						nb_meas_per_KF[k]++;
					}
				}
				
			}
			//std::cout<<"nb meas view ["<<k<<"] = "<<vdErrorSquared.size()<<std::endl;
			sigma_tukey[k]=0;
			if(vdErrorSquared.size()>0)
				sigma_tukey[k]=getSigmaSquared(vdErrorSquared);
			if(sigma_tukey[k] < MINSIGMATUKEY)
				sigma_tukey[k] = MINSIGMATUKEY;	
		}
		
		
		//now collect all jacobians for camera and depth features to optimise (for that needs to be matched, have a positive depth which is not outlier)
		//get LUT for depths to optimise
		//std::cout<<"get LUTs"<<std::endl;
		
		int LUTcam_to_opt[nb_kf];
		int LUTcam_to_main[nb_kf];
		int cam_opt_cpt=0;
		for(int k=0;k<nb_kf;k++)
		{
			if(nb_meas_per_KF[k]>4)
			{
				LUTcam_to_opt[k]=cam_opt_cpt;
				LUTcam_to_main[cam_opt_cpt]=k;
				cam_opt_cpt++;
			}
			else
				LUTcam_to_opt[k]=-1;
		}
		
		//get error and jacobians
		//std::cout<<"get error and jacs"<<std::endl;
		std::vector<miniBAjacz> fJacobian;		

		for(int k=0;k<nb_kf;k++)
		{		
			for(int i=0;i<nb_feat_ref;i++)
				if(std::find(ptToUse.begin(),ptToUse.end(),i)!=ptToUse.end())
			{
				//Vector3f mapPointGT=relPoseGT[k]*(toHomogeneous(invDepthGT[i].meterCoord)/invDepthGT[i].invDepth);
				Vector3f mapPointEstim=relPoseEstim[k]*(toHomogeneous(invDepthEstim[i].meterCoord)/invDepthEstim[i].invDepth);
				//if(mapPointGT[2]>0 && mapPointEstim[2]>0)
				{
					//Vector2f x_d=ProjInvDepthPoint(invDepthGT[i],relPoseGT[k]);//current projection
					//Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
						Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels

						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
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
							newJc.index_cam_opt=LUTcam_to_opt[k];
							if(newJc.index_cam_opt!=-1)
							{
							
								MatrixXf de_dp(2,6);
								de_dp(0,0)=-w/D;	de_dp(0,1)=0;		de_dp(0,2)=Nx*w/(D*D);	
								de_dp(1,0)=0;		de_dp(1,1)=-w/D;	de_dp(1,2)=Ny*w/(D*D);	
								
								de_dp(0,3)=rotCoord[1]*Nx/(D*D);	de_dp(0,4)=(-D-rotCoord[0]*Nx)/(D*D);	de_dp(0,5)=rotCoord[1]/D;
								de_dp(1,3)=(D+rotCoord[1]*Ny)/(D*D);	de_dp(1,4)=-rotCoord[0]*Ny/(D*D);	de_dp(1,5)=-rotCoord[1]/D;
								
								newJc.de_dp=de_dp;
							}
							
							fJacobian.push_back(newJc);

						}
					}
				}
			}
	
		}		
		//compute Hessain and update
		int nbCamsToUpdate=nb_kf;
		VectorXf Jtep[nbCamsToUpdate];
		for(int i=0;i<nbCamsToUpdate;i++)
		{
			Jtep[i].resize(6);
			Jtep[i].setZero();
		}
		MatrixXf Hpp[nbCamsToUpdate];
		for(int i=0;i<nbCamsToUpdate;i++)
		{
			Hpp[i].resize(6,6);
			Hpp[i].setZero();
		}
		
		
		for(int i=0;i<fJacobian.size();i++)			
		{
			short &pt_cam_id=fJacobian[i].index_cam_opt;

			if(pt_cam_id!=-1)
			{
				//update Jte
				Jtep[pt_cam_id]+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].proj_error;
				//update Hessian
				Hpp[pt_cam_id]+=fJacobian[i].weight *fJacobian[i].de_dp.transpose()*fJacobian[i].de_dp;	
			}
		}
		
		
		for(int k=0;k<nbCamsToUpdate;k++)
		for(int i=0;i<6;i++)
			Hpp[k](i,i)=(1.+mLMLambda)*Hpp[k](i,i);	

		
		HomogeneousMatrix newRelPoseEstim[nb_kf];
		for(int k=0;k<nb_kf;k++)newRelPoseEstim[k]=relPoseEstim[k];
		  
		for(int k=0;k<nbCamsToUpdate;k++)
		{
			Eigen::FullPivLU<MatrixXf> lu(Hpp[k]);
			VectorXf Dci=-0.2*(lu.inverse()*Jtep[k]);
			if(isnan(Dci[0])){coutRed<<"nanerie cam"<<endlRed;return -1;};
			short id_glob=LUTcam_to_main[k];
			//std::cout<<"Dci["<<k<<"] = "<<Dci<<std::endl;
			//newRelPoseEstim[k]=HomogeneousMatrix(Dci)*relPoseEstim[k];
			//use rotation in compositional and translkation in additional
			newRelPoseEstim[id_glob].set_rotation(HomogeneousMatrix(Dci).get_rotation()*relPoseEstim[id_glob].get_rotation());
			newRelPoseEstim[id_glob].set_translation(HomogeneousMatrix(Dci).get_translation()+relPoseEstim[id_glob].get_translation());
		}
		
		
		for(int k=0;k<nb_kf;k++)
		{		
			float residue_robust=0;
			float valid_wmeas=0;
			float residue_robustAfter=0;
			float valid_wmeasAfter=0;
			for(int i=0;i<nb_feat_ref;i++)
				if(std::find(ptToUse.begin(),ptToUse.end(),i)!=ptToUse.end())
			{
				{
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f x_c2=ProjInvDepthPoint(invDepthEstim[i],newRelPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
						Vector2f error2=(x_d-x_c2);//error in pixels
						Vector2f errorPix=myCamera.m2PixProjJac()*error;//error in pixels
						Vector2f errorPix2=myCamera.m2PixProjJac()*error2;//error in pixels

						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
							residue_robust+=error.transpose()*error;
							valid_wmeas++;
						}
						else
						{
							residue_robust+=sigma_tukey[k];
							valid_wmeas++;
						}
						float TukeyCoef2=squareRootTukey(errorPix2.transpose()*errorPix2,sigma_tukey[k]);
						if(TukeyCoef2>0)
						{
							residue_robustAfter+=error2.transpose()*error2;
							valid_wmeasAfter++;
						}
						else
						{
							residue_robustAfter+=sigma_tukey[k];
							valid_wmeasAfter++;
						}
					}

				}
			}
			if(valid_wmeasAfter>0 && residue_robustAfter/valid_wmeasAfter<residue_robust/valid_wmeas)
				relPoseEstim[k]=newRelPoseEstim[k];
			
		}

	}

	
	for(int k=0;k<nb_kf;k++)
		PoseEstim[k]=relPoseEstim[k]*PoseRef;
}

int optimPointPose()
{
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;	
	
	HomogeneousMatrix relPoseEstim[nb_kf];
	for(int k=0;k<nb_kf;k++)
	{
		relPoseEstim[k]=PoseEstim[k]*PoseRef.inverse();
	}
	
	for(int iter=0;iter<10;iter++)
	{
		//get Tukey factor
		std::cout<<"iter "<<iter<<std::endl;
		//std::cout<<"get Tukey factor iter "<<iter<<std::endl;
		
		
		//now collect all jacobians for camera and depth features to optimise (for that needs to be matched, have a positive depth which is not outlier)
		//get LUT for depths to optimise
		//std::cout<<"get LUTs"<<std::endl;
		int LUT_to_opt[nb_feat_ref];
		int LUT_to_main[nb_feat_ref];
		
		int opt_cpt=0;
		for(int i=0;i<nb_feat_ref;i++)
		{
			LUT_to_opt[i]=-1;
			for(int k=0;k<nb_kf;k++)
			{
				{
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection

						Vector2f error=myCamera.m2PixProjJac()*(x_d-x_c);//error in pixels

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
							LUT_to_opt[i]=opt_cpt;
							LUT_to_main[opt_cpt]=i;
								
							opt_cpt++;
							break;//check next feature
						}
						
					}
				}
			}
		}

		
		//get error and jacobians
		//std::cout<<"get error and jacs"<<std::endl;
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
					//Vector2f x_d=noisyMeasures[k][i];
					Vector2f x_d=AllMatches[i][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[i],relPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
								
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

						
						fJacobian.push_back(newJc);

						
					}
				}
			}
	
		}
		
		int nbPointsToUpdate=opt_cpt;
		VectorXf Jtex(nbPointsToUpdate);Jtex.setZero();
		float Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i]=0;
		
		for(int i=0;i<fJacobian.size();i++)			
		{
			short &pt_opt_id=fJacobian[i].index_pt_opt;
			if(pt_opt_id!=-1)
			{
				//update Jte
				Jtex[pt_opt_id]+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].proj_error;
				//update Hessian
				Hxx[pt_opt_id]+=fJacobian[i].weight *fJacobian[i].de_dz.transpose()*fJacobian[i].de_dz;
			}
		}
		
		
		for(int i=0;i<nbPointsToUpdate;i++)
			Hxx[i]=(1.+mLMLambda)*Hxx[i];	
		
		

		//update point coord and depths
		for(int id_opt=0;id_opt<nbPointsToUpdate;id_opt++)
		{
			short id_glob=LUT_to_main[id_opt];
			float Dz=0.2*-Jtex[id_opt]/Hxx[id_opt];
			float newIz=invDepthEstim[id_glob].invDepth+Dz;
			PointInvDepth newInvDepthEstim=invDepthEstim[id_glob];
			newInvDepthEstim.invDepth=newIz;
			if(!isnan(newIz) && !isinf(newIz))
			{
				float residue=0;
				float residueAfter=0;
				for(int k=0;k<nb_kf;k++)
				{	
					Vector2f x_d=AllMatches[id_glob][k];
					if(x_d[0]!=-1)
					{
						Vector2f x_c=ProjInvDepthPoint(invDepthEstim[id_glob],relPoseEstim[k]);//current projection
						Vector2f x_c2=ProjInvDepthPoint(newInvDepthEstim,relPoseEstim[k]);//current projection
						Vector2f error=(x_d-x_c);//error in pixels
						Vector2f error2=(x_d-x_c2);//error in pixels
						residue+=error.transpose()*error;
						residueAfter+=error2.transpose()*error2;
					}	
				}	
				
				if(residueAfter<residue)
					invDepthEstim[id_glob].invDepth=newIz;
			}
		}
			
	}
}
