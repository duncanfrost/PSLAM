//load virtual sequence
//do tracking of fast feature to get stereo pair
//Init map from stereo match

#define AMOVERBOSE 1

#include <fstream>
#include "../src/Primitives/Camera.h" 
#include "../src/Visualisation/VisualisationModule.h"
#include "../src/ImageSource/VideoSourceLiveCV.h"
#include "../src/MapEngines/MapOptimiserEssential.h"
#include "../src/TrackEngines/MapTracker.h"
#include "../src/Primitives/obMap.h"

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
SimuWindow *mSimuWindow;

//camera object, calibration ...
Camera myCamera;
//estimated map
obMap Map_Estim;
//obMapRun Map_Estim;
//map tracker
MapTracker *mapTracker;
	

int id_current_frame=0;
HomogeneousMatrix current_estimated_pose;


int main(int argc, char** argv)
{
	InitProcAndGPU();

	std::cout<<"################################################################# "<<std::endl;
	std::cout<<"################    SLAM    SIMU           ###################### "<<std::endl;
	std::cout<<"################################################################# "<<std::endl;
	//get camera calibration
	myCamera.Init(640,480,CamPlaystationEye);
	Map_Estim.InitCam(&myCamera);
	
	//create mapTracker
	mapTracker=new MapTracker(&myCamera,&Map_Estim);	
	
	//create a window to visualize estimated map
	VisuEngine= new VisualisationModule(&Idle);

	VisuEngine->addWindowSimu("Simu",WIMU_WORLD_DIR,"city.obj",640,480,&myCamera,mSimuWindow);
	VisuEngine->setOnDraw(&addDrawFunction);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	
	VisuEngine->addWindowMap("Map",640,480,&myCamera,mMapWindow);
	VisuEngine->setOnKeyPress(&processNormalKeysPlus);
	mMapWindow->addMap(&Map_Estim);
	
	//add current estimated cam position in map viewer
	mMapWindow->addCamera(&current_estimated_pose,Vector3f(1.,0.,0.));
	
	//HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
	HomogeneousMatrix moveCam(0.,-3.,3.,-1.0,0,0);
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

void Idle(void) 
{
	//grab new image
	cv::Mat &current_img_BW=mSimuWindow->getCurrentImageBW();

	//process it
	if(!pause_process)
	{
		//mapTracker->TrackFrame(current_img_BW);
		mapTracker->TrackFrame(current_img_BW,&mSimuWindow->getCurrentImage(),true);
		current_estimated_pose=mapTracker->getPose();//update current camera pose for map viewer
	
		mMapWindow->showClosestKF(mapTracker->getIdRelKF());
		mMapWindow->showActiveKF(mapTracker->getClosestKFs());
		
	}
	
	
	VisuEngine->drawWindows();
	id_current_frame++;
}
void addDrawFunction(void) 
{	
	//std::cout<<"Display"<<std::endl;
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

void processNormalKeysPlus(unsigned char key, int x, int y)
{
	switch(key) {
		case 'm':
			Map_Estim.saveToFile("map.dat");
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

