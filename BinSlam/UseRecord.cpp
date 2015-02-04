//same as live but myVideoSource is just set to use sequence instead of live capture

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceSeq.h"
#include "../src/TrackEngines/MapTracker.h"
#include "../src/Primitives/obMap.h"
#include "../src/MapEngines/BundleAdjuster.h"

#include <Eigen/Core>
using namespace Eigen;

void Idle(void) ;
void processNormalKeysPlus(unsigned char key, int x, int y);//process key on top of existing one from mapViewer
void addDrawFunction(void) ;
void onClickShowTexture(void) ;
void onClickShowMatching(void) ;
void onClickShowLocFeatures(void) ;

//Visualization
VisualisationModule *VisuEngine;
//map viewer from VisuEngine
MapWindow *mMapWindow;
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
	

int id_current_frame=0;
HomogeneousMatrix current_estimated_pose;

PlottingWindow *PlotWindow;

int main(int argc, char** argv)
{

    #ifdef USE_OMP_C
        coutGreen << "We are using openMP" << endlGreen;
    #endif


	InitProcAndGPU();

	std::cout<<"################################################################# "<<std::endl;
	std::cout<<"################    SLAM    NEW            ###################### "<<std::endl;
	std::cout<<"################################################################# "<<std::endl;
	//get camera calibration
	//myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	myVideoSource=new VideoSourceSeq("/home/amaury/Program/Projects/files/record/%06d.png",CamPlaystationEye,0);
	//myVideoSource->setIncrement(-1);

	myVideoSource->initCam();
	
	
	
	myCamera=*myVideoSource->getPointerCamera();
	Map_Estim.InitCam(&myCamera);
	
	//create mapTracker
	mapTracker=new MapTracker(&myCamera,&Map_Estim);	
	
	//create a window to visualize estimated map
	VisuEngine= new VisualisationModule(&Idle);
	VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	
	//VisuEngine->addWindowMap("Map",myCamera.get_width(),myCamera.get_height(),&myCamera,mMapWindow);
	
	Camera mCamDouble(2*myCamera.get_width(),2*myCamera.get_height(),CamPlaystationEye);
	VisuEngine->addWindowMap("Map",2*myCamera.get_width(),2*myCamera.get_height(),&mCamDouble,mMapWindow);
	
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	VisuEngine->addWindowPlot("Confidence last edge",640,480,PlotWindow,1,200,false);
	PlotWindow->setValuesPerPlot(0,1);//just want one plot, with one value to display
	mMapWindow->addMap(&Map_Estim);
	
	//add current estimated cam position in map viewer
	mMapWindow->addCamera(&current_estimated_pose,Vector3f(1.,0.,0.));
	
	//HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
	HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
	mMapWindow->setCameraPose(moveCam);
	
	std::cout<<"Start"<<std::endl;

	VisuEngine->prepareLoop(argc, argv);
	mMapWindow->addButton("Textures",&onClickShowTexture);
	mMapWindow->addButton("Matching",&onClickShowMatching);
	mMapWindow->addButton("LocalFeat",&onClickShowLocFeatures);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


amoTimer time0;
bool pause_process=false;
bool stopAtEachFrame=false;

void Idle(void) 
{
	if(!pause_process)
	{
		//grab new image
		myVideoSource->grabNewFrame();
		cv::Mat &current_img_BW=*myVideoSource->GetFramePointerBW();

		//process it
		//mapTracker->TrackFrame(current_img_BW);
		mapTracker->TrackFrame(current_img_BW,myVideoSource->GetFramePointer(),true);
		current_estimated_pose=mapTracker->getPose();//update current camera pose for map viewer
		if(stopAtEachFrame)pause_process=true;
		
		//get confidence last edge and plot it
		float conf;
		if(Map_Estim.getKF(Map_Estim.getNbKeyFrames()-1)->getNbNeigbours()>0)
		{
			conf=Map_Estim.getKF(Map_Estim.getNbKeyFrames()-1)->getPtNeigbour(0)->edgeScore;
			PlotWindow->setVal(0,0,conf);
			PlotWindow->incrementTimeLine();
		}
		
		mMapWindow->showClosestKF(mapTracker->getIdRelKF());
		mMapWindow->showActiveKF(mapTracker->getClosestKFs());

	}
	
	
	VisuEngine->drawWindows();
	id_current_frame++;

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
		
	std::vector<int> closeKF=mapTracker->getClosestKFs();
	for(int c=0;c<closeKF.size();c++)
	{
		//std::cout<<"display close kf id = "<<closeKF[c]<<std::endl;
		int kfn=closeKF[c];
		KeyFrame &kfc=*Map_Estim.getKF(kfn);
		std::vector<p_match> matches=kfc.getCurrentMatches();
		
		int nbMatchMapPoint=0;
		float residueMapPoint = 0;
		int nbMatchFeature=0;
		float residueFeatures = 0;
		
		for(int i=0;i<matches.size();i++)
		{

			glColor3f(0,1.,0);
			glBegin(GL_POINTS);	
			glVertex2f(matches[i].u1c,matches[i].v1c);
			glEnd();
			
			//project all matched features
			int id_feat=kfc.indexCandidateFeatureFromVisoId(matches[i].i1p);
			if(id_feat!=-1)
			{
				uptoscaleFeature &feat=*kfc.getPtLocalBestFeatures(id_feat);
				if(feat.matched)
				{
					MapPoint &point=*feat.ptKForigin->getPtMapPoint(feat.idPoint);
					Vector3f mapPointsCam=current_estimated_pose*point.getPosition();
					Vector2f projPix=myCamera.Project(mapPointsCam);//current projection
					
					glColor3f(1,0.,0);
					glBegin(GL_LINES);	
					glVertex2f(projPix[0],projPix[1]);
					glVertex2f(matches[i].u1c,matches[i].v1c);
					glEnd();
					
					nbMatchMapPoint++;
					residueMapPoint+=sqrt((projPix[0]-matches[i].u1c)*(projPix[0]-matches[i].u1c)+(projPix[1]-matches[i].v1c)*(projPix[1]-matches[i].v1c));
				}
				else
				{
					Vector3f mapPointsCam=current_estimated_pose*kfc.getPose().inverse()*feat.getLocalCoordinates();
					Vector2f projPix=myCamera.Project(mapPointsCam);//current projection
					
					glColor3f(1,0.6,0);
					glBegin(GL_LINES);	
					glVertex2f(projPix[0],projPix[1]);
					glVertex2f(matches[i].u1c,matches[i].v1c);
					glEnd();	
					
					nbMatchFeature++;
					residueFeatures+=sqrt((projPix[0]-matches[i].u1c)*(projPix[0]-matches[i].u1c)+(projPix[1]-matches[i].v1c)*(projPix[1]-matches[i].v1c));
				}
			}
			
		}
	}
	glColor3f(1,1,1);
	unset2DGLProjection();

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
void updateScaleLastKeyFrame2()
{
	//get list of KF to optimise (everything but first)
	int width_inner=50;
		
	if(width_inner>Map_Estim.getNbKeyFrames()-1)width_inner=Map_Estim.getNbKeyFrames()-1;
	std::vector<int> innerWind;	
	for(int i=Map_Estim.getNbKeyFrames()-width_inner;i<Map_Estim.getNbKeyFrames();i++)
		innerWind.push_back(i);	
	
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
void doClassicBA()
{
	//get list of KF to optimise (everything but first)
	int width_inner=20;
	
	if(width_inner>Map_Estim.getNbKeyFrames()-1)width_inner=Map_Estim.getNbKeyFrames()-1;
	std::vector<int> optimKF;	
	for(int i=Map_Estim.getNbKeyFrames()-width_inner;i<Map_Estim.getNbKeyFrames();i++)
		optimKF.push_back(i);
	
	BundleAdjuster mBA(&Map_Estim);
	mBA.optimiseInnerWindow(optimKF,1,true);
	
}
void doClassicMinFeatureVariance()
{
       std::cout<<"doClassicMinFeatureVariance"<<std::endl;
	//get list of KF to optimise (everything but first)
	int width_inner=20;
	
	if(width_inner>Map_Estim.getNbKeyFrames()-1)width_inner=Map_Estim.getNbKeyFrames()-1;
	std::vector<int> optimKF;	
	for(int i=Map_Estim.getNbKeyFrames()-width_inner;i<Map_Estim.getNbKeyFrames();i++)
		optimKF.push_back(i);
	
	MapOptimiser mBA(&Map_Estim);
	//mBA.optimiseInnerWindow(optimKF);
	mBA.optimiseInnerWindowRobust(optimKF);
}
void doNewMapOptim()
{
        std::cout<<"optimiseInnerWindowRobust"<<std::endl;
	std::vector<int> optimKF;
	for(int i=0;i<Map_Estim.getNbKeyFrames()-1;i++)optimKF.push_back(i);
		
	int nb_iter=1;
	
	MapOptimiserEssential mBA(&Map_Estim);
	mBA.optimiseInnerWindowRobust(optimKF,nb_iter);
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

void processNormalKeysPlus(unsigned char key, int x, int y)
{
	switch(key) {

		case 'c':
			//try change scale: x2 last KF
			doClassicMinFeatureVariance();
			break;
		case 'b'://do bundle adjust on last KF
			doClassicBA();
			break;
		case 'w'://do bundle adjust on last KF
			doNewMapOptim(); 
			break;

		case 'm':
			Map_Estim.saveToFile("map.dat");
			break;

		case 's':
			//try change scale: x2 last KF
			updateScaleLastKeyFrame();
			mapTracker->updateRelPoseAfterMapOptim();
			current_estimated_pose=mapTracker->getPose();
			break;
		case 'd':
			//try change scale: x2 last KF
			updateScaleLastKeyFrame2();
			mapTracker->updateRelPoseAfterMapOptim();
			current_estimated_pose=mapTracker->getPose();
			break;
		case ' ':
			pause_process=!pause_process;
			break;
		case 'r':
			removeUnusedPoints();
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

