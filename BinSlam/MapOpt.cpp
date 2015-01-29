//load estimated map and optimise it

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceSeq.cpp"
#include "../src/TrackEngines/MapTracker.h"
#include "../src/Primitives/obMap.h"
#include "../src/MapEngines/BundleAdjuster.h"
#include "../src/MapEngines/MapOptimiserEssential.h"
#include "../src/MapEngines/PoseGraphOpt.h"

#include <Eigen/Core>
using namespace Eigen;

void Idle(void) ;
void processNormalKeysPlus(unsigned char key, int x, int y);//process key on top of existing one from mapViewer
void addDrawFunction(void) ;
void addDrawFunctionLoop(void) ;
void onClickShowTexture(void) ;
void onClickShowMatching(void) ;
void onClickShowLocFeatures(void) ;


//Visualization
VisualisationModule *VisuEngine;
//map viewer from VisuEngine
MapWindow *mMapWindow;
MotherWindow *InKfWindow;
MotherWindow *MatchWindow;
cv::Mat imageKfDisp;
//image acquisition
//VideoSourceLiveCV *myVideoSource;
VideoSourceSeq *myVideoSource;
//camera object, calibration ...
Camera myCamera;
//estimated map
obMap Map_Estim;
//obMapRun Map_Estim;
//map tracker
MapTracker *mapTracker;
	
//for test loop
cv::Mat ImgDisplayMatch;

int id_current_frame=0;
HomogeneousMatrix current_estimated_pose;

int main(int argc, char** argv)
{
	Vector2f vtest;vtest[0]=2;vtest[1]=2;
	std::cout<<vtest.squaredNorm()<<std::endl;
	InitProcAndGPU();

	std::cout<<"################################################################# "<<std::endl;
	std::cout<<"################    SLAM    NEW            ###################### "<<std::endl;
	std::cout<<"################################################################# "<<std::endl;
	//get camera calibration
	//myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	//myVideoSource=new VideoSourceSeq("/home/amaury/Program/Projects/files/KITTI/%06d.png",KittiCam,0);
	//myVideoSource->initCam();
	//myCamera=*myVideoSource->getPointerCamera();
	myCamera.Init(1241,376,KittiCam);
	//myCamera.Init(640,480,CamPlaystationEye);
	Map_Estim.InitCam(&myCamera);
	
	Map_Estim.loadFromFile("map.dat");
	std::cout<<"nb points = "<<Map_Estim.getNbMapPoints()<<std::endl;
	std::cout<<"nb used points = "<<Map_Estim.getNbUsedMapPoints()<<std::endl;
	
	Map_Estim.getKF(0)->getImg_p(0).copyTo(imageKfDisp);
	
	//move kf 0 to Id
	HomogeneousMatrix pose0=Map_Estim.getKF(0)->getPose();
	for(int i=0;i<Map_Estim.getNbKeyFrames();i++)
	{
		Map_Estim.getKF(i)->setPose(pose0.inverse()*Map_Estim.getKF(i)->getPose());
		for(int j=0;j<Map_Estim.getKF(i)->getNbMapPoint();j++)
		{
			MapPoint &point= *Map_Estim.getKF(i)->getPtMapPoint(j);
			point.updatePosition(pose0*point.getPosition());
		}
	}
	
	ImgDisplayMatch.create(Map_Estim.getKF(0)->getImg_p(lvl_Viso).size().height, Map_Estim.getKF(0)->getImg_p(lvl_Viso).size().width*2, CV_8UC3);

	
	//HomogeneousMatrix errorIntroduced(0.0,0.1,0,0,0,0);
	//Map_Estim.getKF(0)->setPose(errorIntroduced*Map_Estim.getKF(0)->getPose());
	
	//create mapTracker
	mapTracker=new MapTracker(&myCamera,&Map_Estim);	
	
	//create a window to visualize estimated map
	VisuEngine= new VisualisationModule(&Idle);
	//VisuEngine->addWindowEmpty("View KF ",myCamera.get_width(),myCamera.get_height(),&addDrawFunction,InKfWindow,&processNormalKeysPlus);
	VisuEngine->addWindowImage("SurfMatching",&ImgDisplayMatch,MatchWindow);
	VisuEngine->setOnDraw(&addDrawFunctionLoop);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	
	VisuEngine->addWindowImage("View KF ",&imageKfDisp,InKfWindow);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	
	VisuEngine->addWindowMap("Map",myCamera.get_width(),myCamera.get_height(),&myCamera,mMapWindow);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	mMapWindow->addMap(&Map_Estim);
	
	
	//add current estimated cam position in map viewer
	mMapWindow->addCamera(&current_estimated_pose,Vector3f(1.,0.,0.));
	
	//HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
	//HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
	//mMapWindow->setCameraPose(moveCam);
	HomogeneousMatrix moveCam(-0.582611, -1.76823,  1.67988, -1.48444, 0.122459, 0.178899);
	mMapWindow->setCameraPose(moveCam);
	
	std::cout<<"##################################################################"<<std::endl;
	std::cout<<"Map Viewer commands:"<<std::endl;
	std::cout<<"press:\t -\'d\' to disturb keyframe pose"<<std::endl;
	std::cout<<"\t -\'b\' for standard BA"<<std::endl;
	std::cout<<"\t -\'c\' for min variance of points"<<std::endl;
	std::cout<<"\t -\'q\' for essential matrix + scale optim"<<std::endl;
	std::cout<<"\t -\'w\' for pose graph optim"<<std::endl;
	std::cout<<"Loop closure commands:"<<std::endl;
	std::cout<<"\t -\'l\' for do matching between kfs"<<std::endl;
	std::cout<<"\t -\'p\' to create edge between newly matched kfs"<<std::endl;
	std::cout<<"Keyframe Viewer commands:"<<std::endl;
	std::cout<<"press:\t -\'k\' to see next keyframe"<<std::endl;
	std::cout<<"##################################################################"<<std::endl;

	VisuEngine->prepareLoop(argc, argv);
	mMapWindow->addButton("Textures",&onClickShowTexture);
	mMapWindow->addButton("Matching",&onClickShowMatching);
	mMapWindow->addButton("LocalFeat",&onClickShowLocFeatures);
	
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


amoTimer time0;
bool pause_process=false;
bool stopAtEachFrame=true;

void Idle(void) 
{
	
	VisuEngine->drawWindows();
}


void updateScaleLastKeyFrame()
{
	//get list of KF to optimise (everything but first)
	std::vector<int> innerWind;
	for(int i=1;i<Map_Estim.getNbKeyFrames();i++)innerWind.push_back(i);
	
	//create list of new scale constraint:
	std::vector<scale_constraint> mScaleConstraints;
	
	//only constraint to test here is scale of last kf must be 2 times bigger
	scale_constraint scale_constraint1;
	scale_constraint1.kf_id=Map_Estim.getNbKeyFrames()-1;
	scale_constraint1.rescale=2;
	mScaleConstraints.push_back(scale_constraint1);
	
	//do optim
	MapOptimiser mBA(&Map_Estim);
	mBA.optimiseScale(innerWind,mScaleConstraints);
	
	
	
}
void doClassicPGO()//'w'
{
	std::vector<int> optimKF;	
	for(int i=1;i<Map_Estim.getNbKeyFrames();i++)optimKF.push_back(i);
	//for(int i=0;i<Map_Estim.getNbKeyFrames()-1;i++)optimKF.push_back(i);
	
	//Map_Estim.poseGraphInnerWindow(optimKF,1);	
	PoseGraphOptimiser mPGO(&Map_Estim);
	amoTimer timer;
	timer.start();
	mPGO.optimiseInnerWindow(optimKF,10);	
	timer.stop("PGO");
	
}
void doMapOptimEssential()//'q'
{
	
	std::vector<int> optimKF;
	for(int i=1;i<Map_Estim.getNbKeyFrames();i++)optimKF.push_back(i);
	//for(int i=0;i<Map_Estim.getNbKeyFrames()-1;i++)optimKF.push_back(i);
	int nb_iter=1;
	
	MapOptimiserEssential mBA(&Map_Estim);
	//mBA.optimiseInnerWindow(optimKF,nb_iter);
	mBA.optimiseInnerWindowRobust(optimKF,nb_iter);

}

void doClassicBA()//'b'
{
	//get list of KF to optimise (everything but first)
	std::vector<int> optimKF;	
	//for(int i=1;i<Map_Estim.getNbKeyFrames();i++)optimKF.push_back(i);
	//for(int i=0;i<Map_Estim.getNbKeyFrames()-1;i++)optimKF.push_back(i);
	optimKF.push_back(0);
	
	BundleAdjuster mBA(&Map_Estim);
	mBA.optimiseInnerWindow(optimKF,1,true);	
	
}

void doClassicMinFeatureVariance()//'c'
{
	std::vector<int> optimKF;	
	for(int i=1;i<Map_Estim.getNbKeyFrames();i++)optimKF.push_back(i);
	//for(int i=0;i<Map_Estim.getNbKeyFrames()-1;i++)optimKF.push_back(i);
	
	MapOptimiser mBA(&Map_Estim);
	//mBA.optimiseInnerWindow(optimKF);
	amoTimer timer;
	timer.start();
	mBA.optimiseInnerWindow2(optimKF);
	//mBA.optimiseInnerWindowRobust(optimKF);
	timer.stop("MinVar");
	//mBA.optimiseInnerWindowRobust2(optimKF);
	std::cout<<"nb points = "<<Map_Estim.getNbMapPoints()<<std::endl;
	std::cout<<"nb used points = "<<Map_Estim.getNbUsedMapPoints()<<std::endl;
}


void removeUnusedPoints()
{
	for(int i=0;i<Map_Estim.getNbKeyFrames();i++)
	{
		std::cout<<"remove unused points kf["<<i<<"]"<<std::endl;
		Map_Estim.removeUnusedPoints(i);
	}
	std::cout<<"nb points = "<<Map_Estim.getNbMapPoints()<<std::endl;
	std::cout<<"nb used points = "<<Map_Estim.getNbUsedMapPoints()<<std::endl;
}
#include "../src/TrackEngines/RobustMatching.h"
#include <vector_types.h>
std::vector<Vector2f> scene_corners_disp(4);
std::vector<Vector2f> scene_corners_disp_refined(4);
std::vector<p_match> ORBsMatches;
std::vector<p_match> refinedMatches;
int kfCheck_prev=0;
int kfCheckWithCurrent_prev=0;
//int kfCheck=2;
//int kfCheckWithCurrent=30;
int kfCheck=34;
int kfCheckWithCurrent=458;

void checkLoopClosure()
{
	refinedMatches.clear();
	/*int idCheck=5;
	int depth_min=20;
	int depth_max=60;
	Map_Estim.checkForSmallLoop(idCheck,depth_min,depth_max);*/
	/*int idCheck=0;
	int depth_min=7;
	int depth_max=15;
	Map_Estim.checkForSmallLoop(idCheck,depth_min,depth_max);*/
	
	char newTitle[200];
	sprintf(newTitle,"New Matching between KF %d and %d",kfCheck, kfCheckWithCurrent);
	MatchWindow->setTitle(newTitle);
	std::cout<<"check "<<kfCheck << " and "<< kfCheckWithCurrent<<std::endl;
	
	KeyFrame &KF1=*Map_Estim.getKF(kfCheck);
	KeyFrame &KF2=*Map_Estim.getKF(kfCheckWithCurrent);

	cv::Mat &Img1=KF1.getImg_p(lvl_Viso);
	cv::Mat &Img2=KF2.getImg_p(lvl_Viso);
	
	cv::Mat Img1Col;cv::cvtColor(Img1,Img1Col,CV_GRAY2RGB);
	cv::Mat Img2Col;cv::cvtColor(Img2,Img2Col,CV_GRAY2RGB);
	
	for(int y=0;y<Img1.size().height;y++)
		for(int x=0;x<Img1.size().width;x++)
			ImgDisplayMatch.at<uchar3>(y, x)=Img1Col.at<uchar3>(y, x);	
	//display ref and matches
	for(int y=0;y<Img1.size().height;y++)
		for(int x=0;x<Img1.size().width;x++)
			ImgDisplayMatch.at<uchar3>(y, x+Img2.size().width)=Img2Col.at<uchar3>(y, x);
	
	
	Matrix3f Homography_viso;
	ORBsMatches=checkForLoop(Img1,Img2,Homography_viso);
	Matrix3f Homography=Homography_viso;
	float div_lvl=ScaleLevel(lvl_Viso);
	for(int i=0;i<2;i++)Homography(i,2)=Homography_viso(i,2)*div_lvl;
	for(int i=0;i<2;i++)Homography(2,i)=Homography_viso(2,i)/div_lvl;
	//for(int i=0;i<2;i++)Homography(i,2)=Homography_viso(i,2)/div_lvl;
	//for(int i=0;i<2;i++)Homography(2,i)=Homography_viso(2,i)*div_lvl;
	
	if(ORBsMatches.size()>4)
	{
		std::cout<<"\tORBsMatches.size() = "<<ORBsMatches.size()<<std::endl;
	
			
		std::vector<Vector2f> obj_corners(4);
		obj_corners[0] = Vector2f(0,0);
		obj_corners[1] = Vector2f( Img2.size().width, 0 );
		obj_corners[2] = Vector2f( Img2.size().width, Img2.size().height );
		obj_corners[3] = Vector2f( 0,  Img2.size().height);
		
		for(int i=0;i<4;i++)
		{
			float x=obj_corners[i][0];
			float y=obj_corners[i][1];
			
			//double denom=Homography(2,0)*x+Homography(2,1)*y+Homography(2,2);
			//scene_corners_disp[i][0]=(Homography(0,0)*x+Homography(0,1)*y+Homography(0,2))/denom;
			//scene_corners_disp[i][1]=(Homography(1,0)*x+Homography(1,1)*y+Homography(1,2))/denom;
			double denom=Homography_viso(2,0)*x+Homography_viso(2,1)*y+Homography_viso(2,2);
			scene_corners_disp[i][0]=(Homography_viso(0,0)*x+Homography_viso(0,1)*y+Homography_viso(0,2))/denom;
			scene_corners_disp[i][1]=(Homography_viso(1,0)*x+Homography_viso(1,1)*y+Homography_viso(1,2))/denom;
		}	
		
		KF1.setHomography(Homography);
		KF1.computeOverlap();
		float overlap=KF1.getOverlapWithLastFrame();
		std::cout<<"\tHomography found with overlap = "<<overlap<<std::endl;
		
		int nb_matchs=KF1.useNewFrame(KF2.getImg_p(),&myCamera);
		refinedMatches=KF1.getCurrentMatches();
		std::cout<<"\trefinedMatches = "<<refinedMatches.size()<<std::endl;
		
		for(int i=0;i<refinedMatches.size();i++)
		{
			refinedMatches[i].u1p=LevelNPos(refinedMatches[i].u1p,lvl_Viso);
			refinedMatches[i].v1p=LevelNPos(refinedMatches[i].v1p,lvl_Viso);
			refinedMatches[i].u1c=LevelNPos(refinedMatches[i].u1c,lvl_Viso);
			refinedMatches[i].v1c=LevelNPos(refinedMatches[i].v1c,lvl_Viso);
		}		
		
		Homography=KF1.getHomography();
		for(int i=0;i<2;i++)Homography_viso(i,2)=Homography(i,2)/div_lvl;
		for(int i=0;i<2;i++)Homography_viso(2,i)=Homography(2,i)*div_lvl;
		//for(int i=0;i<2;i++)Homography_viso(i,2)=Homography(i,2)*div_lvl;
		//for(int i=0;i<2;i++)Homography_viso(2,i)=Homography(2,i)/div_lvl;
		KF1.computeOverlap();
		float newOverlap=KF1.getOverlapWithLastFrame();		
		std::cout<<"\tHomography refined with overlap = "<<overlap<<std::endl;
		
		
		for(int i=0;i<4;i++)
		{
			float x=obj_corners[i][0];
			float y=obj_corners[i][1];
			
			//double denom=Homography(2,0)*x+Homography(2,1)*y+Homography(2,2);
			//scene_corners_disp_refined[i][0]=(Homography(0,0)*x+Homography(0,1)*y+Homography(0,2))/denom;
			//scene_corners_disp_refined[i][1]=(Homography(1,0)*x+Homography(1,1)*y+Homography(1,2))/denom;
			double denom=Homography_viso(2,0)*x+Homography_viso(2,1)*y+Homography_viso(2,2);
			scene_corners_disp_refined[i][0]=(Homography_viso(0,0)*x+Homography_viso(0,1)*y+Homography_viso(0,2))/denom;
			scene_corners_disp_refined[i][1]=(Homography_viso(1,0)*x+Homography_viso(1,1)*y+Homography_viso(1,2))/denom;
		}	
		
	}
	else
	{	
		for(int i=0;i<4;i++)scene_corners_disp_refined[i]=Vector2f(-1,-1);
		for(int i=0;i<4;i++)scene_corners_disp[i]=Vector2f(-1,-1);
	}
	kfCheck_prev=kfCheck;
	kfCheckWithCurrent_prev=kfCheckWithCurrent;
	
	//if(kfCheckWithCurrent==Map_Estim.getNbKeyFrames()-1)kfCheck++;
	//kfCheckWithCurrent=(kfCheckWithCurrent+1)%Map_Estim.getNbKeyFrames();
	kfCheck++;
	
	if(kfCheck>40)
	{
		kfCheck=34;
		kfCheckWithCurrent++;
	}
}

void closeLoop()
{
	Map_Estim.createNewEdge(kfCheck_prev,kfCheckWithCurrent_prev);
	
}

void addDrawFunctionLoop(void) 
{	
	
	//get max response feature
	cv::Mat &ImgFirst=Map_Estim.getKF(0)->getImg_p(lvl_Viso);
	
	set2DGLProjection();
	glPointSize(3.0);
	//glColor3f(0,1,0);
	cv::RNG& rng=cv::theRNG();
	glLineWidth(2.);
	for(int i=0;i<ORBsMatches.size();i++)
	{
		glColor3f(rng(256)/256.,rng(256)/256.,rng(256)/256.);
		glBegin(GL_LINES);
		glVertex2f(ORBsMatches[i].u1p,ORBsMatches[i].v1p);
		glVertex2f(ORBsMatches[i].u1c+ImgFirst.size().width,ORBsMatches[i].v1c);
		glEnd();
	}
	
	glLineWidth(1.);
	for(int i=0;i<refinedMatches.size();i++)
	{
		glColor3f(0,1,0);
		glBegin(GL_LINES);
		glVertex2f(refinedMatches[i].u1p,refinedMatches[i].v1p);
		glVertex2f(refinedMatches[i].u1c+ImgFirst.size().width,refinedMatches[i].v1c);
		glEnd();
	}
	
	glLineWidth(4.);
	//green = homog from orbs
	glColor3f(0,1.,0);
	for(int i=0;i<4;i++)
	{
		glBegin(GL_LINES);
		glVertex2f(scene_corners_disp[i][0]+ImgFirst.size().width,scene_corners_disp[i][1]);
		glVertex2f(scene_corners_disp[(i+1)%4][0]+ImgFirst.size().width,scene_corners_disp[(i+1)%4][1]);
		glEnd();
	}
	
	//blue = homog from viso
	glLineWidth(4.);
	glColor3f(0,0.,1);
	for(int i=0;i<4;i++)
	{
		glBegin(GL_LINES);
		glVertex2f(scene_corners_disp_refined[i][0]+ImgFirst.size().width,scene_corners_disp_refined[i][1]);
		glVertex2f(scene_corners_disp_refined[(i+1)%4][0]+ImgFirst.size().width,scene_corners_disp_refined[(i+1)%4][1]);
		glEnd();
	}
	
	glColor3f(1,1,1);
	unset2DGLProjection();
}
void disturb()
{
	VectorXf poseNoise(6);
	poseNoise.setZero();
	srand (time(NULL));
	for(int j=1;j<Map_Estim.getNbKeyFrames();j++)
	{
	for(int i=0;i<3;i++)poseNoise[i]+=0.02*((double)rand()/(double)RAND_MAX-0.5);
	
	//poseNoise[3]+=0.02*((double)rand()/(double)RAND_MAX);
	for(int i=0;i<3;i++)poseNoise[i+3]+=0.02*((double)rand()/(double)RAND_MAX-0.5);
	for(int i=0;i<3;i++)poseNoise[i+3]+=0.05*((double)rand()/(double)RAND_MAX-0.5);
	
	
	Map_Estim.getKF(j)->setPose(HomogeneousMatrix(poseNoise)* Map_Estim.getKF(j)->getPose());
	//Map_Estim.getKF(1)->setPose(HomogeneousMatrix(poseNoise)* Map_Estim.getKF(1)->getPose());
	//Map_Estim.getKF(0)->setPose(HomogeneousMatrix(poseNoise)* Map_Estim.getKF(0)->getPose());
	}

}
int kf_disp=0;
int disp_best_local=0;

void processNormalKeysPlus(unsigned char key, int x, int y)
{
	switch(key) {


		case 'k':
			char fileName[200];
			if(disp_best_local==0)
			{
				disp_best_local=1;
				sprintf(fileName,"View from best pair kf %d, ",kf_disp);
				//imageKfDisp=&Map_Estim.getKF(kf_disp)->getBestImgPair();
				Map_Estim.getKF(kf_disp)->getBestImgPair().copyTo(imageKfDisp);

			}
			else
			{
				disp_best_local=0;
				kf_disp= (kf_disp+1)% Map_Estim.getNbKeyFrames();
				sprintf(fileName,"View from kf %d, ",kf_disp);
				//imageKfDisp=&Map_Estim.getKF(kf_disp)->getImg_p(0);
				Map_Estim.getKF(kf_disp)->getImg_p(0).copyTo(imageKfDisp);
			}
			InKfWindow->setTitle(fileName);
			break;
		case 'l':
			//check matches
			checkLoopClosure();
			break;
		case 'p':
			closeLoop();//close loop between two checked KF
			break;
		case 'd'://esc
			disturb();
			break;
		case 'b'://esc
			doClassicBA();
			break;
		case 'c'://esc
			doClassicMinFeatureVariance();
			break;
		case 'q'://esc
			doMapOptimEssential();
			break;
		case 'w'://esc
			doClassicPGO();
			break;
		case 's':
			//try change scale: x2 last KF
			updateScaleLastKeyFrame();
			break;
		case 'r':
			removeUnusedPoints();
			break;
		case ' ':
			pause_process=!pause_process;
			break;
		case 27://esc
			//Map_Estim.stopThread();
			delete mapTracker;
			delete VisuEngine;
			exit(0);
			break;
	}


	glutPostRedisplay();
}
//draw map viewed from kf_disp if disp_best_local ==0
//else draw map from best ministereo pair of kf_disp
void addDrawFunction(void) 
{	
	KeyFrame &kfd=*Map_Estim.getKF(kf_disp);
  
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glEnable(GL_POINT_SMOOTH);
	
	set2DGLProjection();
	glColor3f(0,1,0);
	
	if(disp_best_local==0)
	{
		//show all local features (their 3D pos is along measure=> no error to display from here)
		//feature position in green
		//projection of corresponding 3D point matked as red dot
		//error marked as red line
		for(int i=0;i<kfd.getNbLocalBestFeatures();i++)
		{
			uptoscaleFeature &cand=*kfd.getPtLocalBestFeatures(i);
			Vector2f posPix=myCamera.ToPixels(cand.posRef);
			glPointSize(3.0);
			glColor3f(0,1.,0);
			glBegin(GL_POINTS);	
			glVertex2f(posPix[0],posPix[1]);
			glEnd();
			
			//if matched to point then display point projection and error
			glColor3f(1,0.,0);
			glLineWidth(2.);
			if(cand.matched)
			{
				MapPoint &point=*cand.ptKForigin->getPtMapPoint(cand.idPoint);
				//project it
				Vector3f coordInKf=kfd.getPose()*point.getPosition();
				Vector2f proj=myCamera.Project(coordInKf);
				
				glPointSize(2.0);
				glBegin(GL_POINTS);	
				glVertex2f(proj[0],proj[1]);
				glEnd();
				
				glBegin(GL_LINES);	
				glVertex2f(posPix[0],posPix[1]);
				glVertex2f(proj[0],proj[1]);
				glEnd();				
			}
		}	
	}
	else
	{
		//show local stereo pair reconstruction error
		HomogeneousMatrix poseBestRel=kfd.getBestRelPose()*kfd.getPose();
		//=> get best matches
		std::vector<p_match> &bestMatches=*kfd.getBestLocalMatches();
		for(int m=0;m<bestMatches.size();m++)
		{
			//draw measure in image
		  	glPointSize(2.0);
			glColor3f(0,1,0);
			glBegin(GL_POINTS);	
			glVertex2f(bestMatches[m].u1c,bestMatches[m].v1c);
			glEnd();
  
		  
		  
			//get corresponding point:
			int id_feat=kfd.indexCandidateFeatureFromVisoId(bestMatches[m].i1p);
			if(id_feat==-1)continue;
			uptoscaleFeature &feat=*kfd.getPtLocalBestFeatures(id_feat);
		
			if(feat.matched)
			{
				glColor3f(1,0,0);
				//if feature matched then need to link measure to matched point 
				MapPoint &point=*feat.ptKForigin->getPtMapPoint(feat.idPoint);
				Vector3f coordInKf=poseBestRel*point.getPosition();
				Vector2f proj=myCamera.Project(coordInKf);
				
				glBegin(GL_POINTS);	
				glVertex2f(proj[0],proj[1]);
				glEnd();
				
				glBegin(GL_LINES);	
				glVertex2f(bestMatches[m].u1c,bestMatches[m].v1c);
				glVertex2f(proj[0],proj[1]);
				glEnd();				
			}
			else
			{
				glColor3f(0,0,1);
				//if not show pos of local feature
				Vector3f coordInKf=kfd.getBestRelPose()*feat.getLocalCoordinates();
				Vector2f proj=myCamera.Project(coordInKf);
				
				glBegin(GL_POINTS);	
				glVertex2f(proj[0],proj[1]);
				glEnd();
				
				glBegin(GL_LINES);	
				glVertex2f(bestMatches[m].u1c,bestMatches[m].v1c);
				glVertex2f(proj[0],proj[1]);
				glEnd();				
			}

			
		}

	}
	glColor3f(1,1,1);
	unset2DGLProjection();
	//glutSwapBuffers();
}

void onClickShowTexture()
{
	std::cout<<"Show textures"<<std::endl;
	mMapWindow->showTextures();
}
void onClickShowMatching()
{
	std::cout<<"Show Matching"<<std::endl;
	mMapWindow->showFeatureConnections();
}
void onClickShowLocFeatures()
{
	std::cout<<"Show Local Features"<<std::endl;
	mMapWindow->showLocalFeature(); 
}

