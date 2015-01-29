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
#include "../src/MapEngines/BundleAdjuster.h"
//#include "../src/Primitives/HomogeneousMatrix.h"


#include "../src/AmoDefines.h"
#include "../src/Primitives/Camera.h"
//#include "ImageProcess.h"
#include "../src/Primitives/obMap.h"
#include "../src/Primitives/HomogeneousMatrix.h"

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
VideoSourceLiveCV *myVideoSource;
//VideoSourceSeq *myVideoSource;
//camera object, calibration ...
Camera myCamera;
//estimated map
obMap Map_Estim;
//obMapRun Map_Estim;
//map tracker
MapTracker *mapTracker;

using namespace std;
	

int id_current_frame=0;
HomogeneousMatrix22 current_estimated_pose;


int main(int argc, char** argv)
{
	InitProcAndGPU();


    //get camera calibration

	myVideoSource=new VideoSourceLiveCV(CamPlaystationEye);
	//myVideoSource=new VideoSourceSeq("/media/Data/Datas/VirtualSeq/new/img_%d.png",CamPlaystationEye,0);
	//myVideoSource->setRecord("/home/amaury/Program/Projects/files/record/%06d.png");
	myVideoSource->initCam();
	myCamera=*myVideoSource->getPointerCamera();
	
	Map_Estim.InitCam(&myCamera);
	

    mapTracker=new MapTracker(&myCamera,&Map_Estim);

    coutRed<< sizeof(Matrix4f) << endlRed;
    coutRed<< sizeof(Matrix4f) << endlRed;

    std::cout << "Address of maptracker: " << mapTracker << std::endl;

    vector<int> vClosest = mapTracker->getClosestKFs();
    cout << sizeof(HomogeneousMatrix22) << endl;


    HomogeneousMatrix22 H;
      cout << sizeof(H) << endl;


    //H.get_HomogMatrix();


//      coutBlue << "Function pointer: " << (void *)Idle << endlBlue;
//    //create a window to visualize estimated map
//    VisuEngine= new VisualisationModule(&Idle);
//    VisuEngine->addWindowImage("Input image",(amoVideoSource*)myVideoSource);
//    VisuEngine->setOnDraw(&addDrawFunction);
//    VisuEngine->setOnKeyPress(&processNormalKeysPlus);

//    VisuEngine->addWindowMap("Map",640,480,&myCamera,mMapWindow);
//    VisuEngine->setOnKeyPress(&processNormalKeysPlus);
//    mMapWindow->addMap(&Map_Estim);
	
//    //add current estimated cam position in map viewer
//    mMapWindow->addCamera(&current_estimated_pose,Vector3f(1.,0.,0.));
	
//    //HomogeneousMatrix moveCam(0.,3.,3.,1.0,0,0);
//    HomogeneousMatrix moveCam(0.,-3.,3.,-1.0,0,0);
//    mMapWindow->setCameraPose(moveCam);

	
//    std::cout<<"##################################################################"<<std::endl;
//    std::cout<<"Commands:"<<std::endl;
//    std::cout<<"press:\t -\'m\' to save map (can be loaded in MapOpt)"<<std::endl;
//    std::cout<<"##################################################################"<<std::endl;

//    VisuEngine->prepareLoop(argc, argv);
//    mMapWindow->addButton("Textures",&onClickShowTexture);
//    mMapWindow->addButton("Matching",&onClickShowMatching);
//    mMapWindow->addButton("LocalFeat",&onClickShowLocFeatures);

//    std::cout << "got here!" << std::endl;

//    coutGreen << mapTracker << endlGreen;
	
//    VisuEngine->startLoop(argc, argv);




	
//	return 0;
}


amoTimer time0;
bool pause_process=false;

void Idle(void) 
{
//   vector<int> vClosest = mapTracker->getClosestKFs();
//   coutRed << mapTracker << endlRed;







//    //grab new image
//    myVideoSource->grabNewFrame();
//    cv::Mat &current_img_BW=*myVideoSource->GetFramePointerBW();
	  
//    //process it
//    if(!pause_process)
//    {
//        mapTracker->TrackFrame(current_img_BW);
//       current_estimated_pose=mapTracker->getPose();//update current camera pose for map viewer
//       coutGreen << current_estimated_pose.get_HomogMatrix() << endlGreen;


	
//        mMapWindow->showClosestKF(mapTracker->getIdRelKF());

//        //vector<int> vClosest = mapTracker->getClosestKFs();
//        mMapWindow->showActiveKF(vClosest);
		
//    }
	
	
    //VisuEngine->drawWindows();
    //id_current_frame++;
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

	//for all active keyframes, show matches from their reference image to current image
	//show the reprojection error of the local 3D features in orange and reproj error of 
	//map points in red
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
			if(id_feat!=-1)//match correspond to a 3D local feature in the keyframe
			{
				uptoscaleFeature &feat=*kfc.getPtLocalBestFeatures(id_feat);
				if(feat.matched)//feature is linked to map point => show error with reprojection of map point
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
				else //not linked to map point=> show reproj error of local feature
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
		case 'r':
			for(int i=0;i<Map_Estim.getNbKeyFrames();i++)
				Map_Estim.removeUnusedPoints(i);
			break;
	
		case 27://esc
			//Map_Estim.stopThread();
			delete mapTracker;
			delete VisuEngine;
			exit(0);
			break;
	}
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

