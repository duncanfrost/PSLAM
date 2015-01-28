//create GT map and create keyframes with corresponding measures
//the noise applied to the measures can be change
//the keyframe poses can be altered to check BA optimisations

#define AMOVERBOSE 1

#include <fstream>
#include <objectr3d/Camera.h>

#include <objectr3d/VisualisationModule.h>
#include <objectr3d/BundleAdjuster.h>
#include <objectr3d/PoseGraphOpt.h>
#include <objectr3d/MapOptimiser.h>
#include <objectr3d/MapOptimiserEssential.h>
#include <objectr3d/obMap.h>

void Idle(void) ;
//load pair of pose and 3D points from real scene
void createMap(obMap &_map);
void optimiseAndOutputForPlotting();


//Visualization
VisualisationModule *VisuEngine;

Camera myCamera;//camera object, calibration ...

//reprojection function and user interface
void keyboard_process(unsigned char key, int x, int y);
void addDrawFunction();
void disturb();

int nb_kf=2;
int nbIter=1;
obMap Map_Estim;

EmptyWindow *InKfWindow;
MapWindow *mMapWindow;
PlottingWindow *PlotWindow;

HomogeneousMatrix PoseGT[100];

//warning: if image noise is big, then algebraic solution for local reconstruction is pretty crap
float image_noise=0.6;

int im_w;
int im_h;

int main(int argc, char** argv)
{
	//Vector2f v;
	//v[0]=2;v[1]=0;
	//std::cout<<v.squaredNorm()<<std::endl;
	im_w=640;
	im_h=480;
	myCamera.Init(im_w,im_h,CamPlaystationEye);
	Map_Estim.InitCam(&myCamera);
	
	//load GT map
	std::cout<<"create map"<<std::endl;
	createMap(Map_Estim);
	std::cout<<"end create map"<<std::endl;
	//disturb();

	
	VisuEngine= new VisualisationModule(&Idle);
	//create a window to visualize reprojection and to grab keyboard event (see keyboard_process)
	VisuEngine->addWindowEmpty("Key frame 0 view",640,480,InKfWindow);
	VisuEngine->setOnDraw(addDrawFunction);
	VisuEngine->setOnKeyPress(keyboard_process);
	//create a window to visualize GT map and estimated one
	VisuEngine->addWindowMap("Map",640,480,&myCamera,mMapWindow);
	VisuEngine->setOnKeyPress(keyboard_process);
	mMapWindow->addMap(&Map_Estim);
	
	VisuEngine->addWindowPlot("PoseError",640,480,PlotWindow,2,200,false);
	PlotWindow->setValuesPerPlot(0,2);//rotation and translation error
	PlotWindow->setValuesPerPlot(1,1);//essential error
	
	
	//HomogeneousMatrix moveCam(0,0.8,0.8,0.7,0,0);
	HomogeneousMatrix moveCam(0,2.1,2.1,2.0,0,0);
	mMapWindow->moveCamera(moveCam);
	
	std::cout<<"##################################################################"<<std::endl;
	std::cout<<"Map Viewer commands:"<<std::endl;
	std::cout<<"press:\t -\'d\' to disturb keyframe pose"<<std::endl;
	std::cout<<"\t -\'b\' for standard BA"<<std::endl;
	std::cout<<"\t -\'c\' for min variance of points"<<std::endl;
	std::cout<<"\t -\'q\' for essential matrix + scale optim"<<std::endl;
	std::cout<<"\t -\'w\' for pose graph optim"<<std::endl;
	std::cout<<"Keyframe Viewer commands:"<<std::endl;
	std::cout<<"press:\t -\'k\' to see next keyframe"<<std::endl;
	std::cout<<"##################################################################"<<std::endl;
	
	VisuEngine->prepareLoop(argc, argv);
	VisuEngine->startLoop(argc, argv);
	
	return 0;
}


void Idle(void) 
{
	//compute pose error and display it
	float error_trans=0;
	float error_rot=0;
	for(int i=0;i<nb_kf;i++)
	{
		HomogeneousMatrix errorPose=PoseGT[i]*Map_Estim.getKF(i)->getPose().inverse();
		VectorXf error_p=errorPose.get_p();
		//std::cout<<"error_p = "<<error_p.transpose()<<std::endl;
		error_trans+=sqrt(error_p.segment(0,3).squaredNorm());
		error_rot+=sqrt(error_p.segment(3,3).squaredNorm());
	}
	
	PlotWindow->setVal(0,0,error_trans);
	PlotWindow->setVal(0,1,error_rot);
	
	float residueEssential=0;
	KeyFrame &KFc=*Map_Estim.getKF(0);
	for(int n=0;n<KFc.getNbNeigbours();n++)
	{
		NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
		KeyFrame &KFn=*Map_Estim.getKF(neigbor.neighboring_kf);
		//compute essential matrix
		HomogeneousMatrix pose_c=KFc.getPose();
		HomogeneousMatrix pose_n=KFn.getPose();
		
		//HomogeneousMatrix relPose=neigbor.relative_poses;
		HomogeneousMatrix relPose=pose_c*pose_n.inverse();
		Vector3f t_nc=relPose.get_translation();
		Matrix3f R_nc=relPose.get_rotation();
		Vector3f nt_nc=t_nc/sqrt(t_nc.squaredNorm());
		//Xc=relPose * Xn;
		//Matrix3f Ecn=GetSkew(t_nc)*relPose.get_rotation();
		Matrix3f Ecn=GetSkew(nt_nc)*relPose.get_rotation();
		
		for(int m=0;m<neigbor.matches.size();m++)
		{
			p_match &mMatch=neigbor.matches[m];
			//measure in current KF KFc
			Vector2f measure_n=myCamera.ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
			//measure in neigbor
			Vector2f measure_c=myCamera.ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
			
			Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
			Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
			
			float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
			//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
			residueEssential+=errorEssential*errorEssential;	
		}
	}
	PlotWindow->setVal(1,0,residueEssential);
	
	PlotWindow->incrementTimeLine();
	
	VisuEngine->drawWindows();
}



void doMapOptimEssential()//'q'
{
	
	std::vector<int> optimKF;
	for(int i=1;i<nb_kf;i++)optimKF.push_back(i);
	//for(int i=0;i<nb_kf-1;i++)optimKF.push_back(i);
	
	MapOptimiserEssential mBA(&Map_Estim);
	mBA.optimiseInnerWindow(optimKF,nbIter);
	//mBA.optimiseInnerWindowRobust(optimKF,nb_iter);
	//mBA.poseGraphRobust(optimKF,nb_iter);

}

void doClassicBA()
{
	//get list of KF to optimise (everything but first)
	std::vector<int> optimKF;	
	for(int i=1;i<nb_kf;i++)optimKF.push_back(i);
	//for(int i=0;i<nb_kf-1;i++)optimKF.push_back(i);
	//optimKF.push_back(0);
	
	BundleAdjuster mBA(&Map_Estim);
	//mBA.optimiseInnerWindow(optimKF,nbIter,true);	
	mBA.optimiseInnerWindow(optimKF,nbIter,false);	
	
}
void doClassicPGO()//'w'
{
	std::vector<int> optimKF;	
	for(int i=1;i<nb_kf;i++)optimKF.push_back(i);
	//for(int i=0;i<nb_kf-1;i++)optimKF.push_back(i);
	
	//Map_Estim.poseGraphInnerWindow(optimKF,1);	
	PoseGraphOptimiser mBA(&Map_Estim);
	//mBA.optimiseInnerWindow(optimKF,nbIter);	
	mBA.optimiseInnerWindow(optimKF,10);	
	
	
}

void doClassicMinFeatureVariance()//c
{
	std::vector<int> optimKF;	
	for(int i=1;i<nb_kf;i++)optimKF.push_back(i);
	//for(int i=0;i<nb_kf-1;i++)optimKF.push_back(i);
	
	MapOptimiser mBA(&Map_Estim);
	mBA.optimiseInnerWindow(optimKF,nbIter);
	//mBA.optimiseInnerWindowRobust(optimKF);
	std::cout<<"nb points = "<<Map_Estim.getNbMapPoints()<<std::endl;
	std::cout<<"nb used points = "<<Map_Estim.getNbUsedMapPoints()<<std::endl;
}

void doClassicMinFeatureVariance2()//c
{
	std::vector<int> optimKF;	
	for(int i=1;i<nb_kf;i++)optimKF.push_back(i);
	//for(int i=0;i<nb_kf-1;i++)optimKF.push_back(i);
	
	MapOptimiser mBA(&Map_Estim);
	mBA.optimiseInnerWindow2(optimKF);
	//mBA.optimiseInnerWindowRobust2(optimKF);
	std::cout<<"nb points = "<<Map_Estim.getNbMapPoints()<<std::endl;
	std::cout<<"nb used points = "<<Map_Estim.getNbUsedMapPoints()<<std::endl;
}

void disturb()
{
	VectorXf poseNoise(6);
	poseNoise.setZero();
	srand (time(NULL));
	for(int j=1;j<nb_kf;j++)
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

void disturb2()
{
	float scale_change=1.5;
	int id_kf_modif=1;
	KeyFrame &KFc=*Map_Estim.getKF(id_kf_modif);
	for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
		KFc.getPtLocalBestFeatures(f)->depthInRef*=scale_change;

	//need to change relative scale in edge and corresponding information matrix
	for(int n=0;n<KFc.getNbNeigbours();n++)
	{
		NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
		Map_Estim.getRelativePoseAndScale(id_kf_modif,neigbor.neighboring_kf,neigbor.relative_scale,neigbor.relative_poses,neigbor.Informationscale,neigbor.InformationMatrixPose);
		
		
		//get edge going in other direction
		KeyFrame &KFn=*Map_Estim.getKF(neigbor.neighboring_kf);
		for(int n2=0;n2<KFn.getNbNeigbours();n2++)
		{
			NeigbourKFNew &neigborinv=*KFn.getPtNeigbour(n2);
			if(neigborinv.neighboring_kf==id_kf_modif)
				Map_Estim.getRelativePoseAndScale(neigbor.neighboring_kf,id_kf_modif,neigborinv.relative_scale,neigborinv.relative_poses,neigborinv.Informationscale,neigborinv.InformationMatrixPose);
			  
		}

	}
	
	HomogeneousMatrix relBestPose= KFc.getBestRelPose();
	relBestPose.set_translation(relBestPose.get_translation()*scale_change);
	KFc.setBestRelPose(relBestPose);
}

void testRelativeConstraint()
{
	MapOptimiser mBA(&Map_Estim);
	float optRelScale;
	HomogeneousMatrix optRelPose;
	float infoScale;
	MatrixXf infoPose;
	
	mBA.getRelativePoseAndScale(0,1,optRelScale,optRelPose,infoScale,infoPose);

	coutRed<<"relPoseInfo"<<endlRed;
	std::cout<<"opt_scale = "<<optRelScale<<std::endl;
	std::cout<<"optRelPose = "<<optRelPose<<std::endl;
	std::cout<<"infoScale = "<<infoScale<<std::endl;
	std::cout<<"infoPose = "<<infoPose<<std::endl;
	
	mBA.getRelativePoseAndScale(1,0,optRelScale,optRelPose,infoScale,infoPose);

	coutRed<<"relPoseInfoInv"<<endlRed;
	std::cout<<"opt_scale = "<<optRelScale<<std::endl;
	std::cout<<"optRelPose = "<<optRelPose<<std::endl;
	std::cout<<"infoScale = "<<infoScale<<std::endl;
	std::cout<<"infoPose = "<<infoPose<<std::endl;
	
	
}

int kf_disp=0;
int disp_best_local=0;

void keyboard_process(unsigned char key, int x, int y)
{
	switch(key) {
		case 'k':
			char fileName[200];
			if(disp_best_local==0)
			{
				disp_best_local=1;
				sprintf(fileName,"View from best pair kf %d, ",kf_disp);

			}
			else
			{
				disp_best_local=0;
				kf_disp= (kf_disp+1)% Map_Estim.getNbKeyFrames();
				sprintf(fileName,"View from kf %d, ",kf_disp);
			}
			InKfWindow->setTitle(fileName);
			break;
		case 27://esc
			exit(0);
			break;
		case 'b'://esc
			doClassicBA();
			break;
		case 'c'://esc
			doClassicMinFeatureVariance();
			//doClassicMinFeatureVariance2();
			break;
		case 'q'://esc
			doMapOptimEssential();
			break;
		case 'd'://esc
			//disturb();
			disturb2();
			break;
		case 'w'://esc
			doClassicPGO();
			break;
		case 'f':
			optimiseAndOutputForPlotting();
			break;
		case 'a':
			testRelativeConstraint();
			break;
			
	}


	glutPostRedisplay();
}



void createMap(obMap &_map)
{
// 	//create loop around point at distance z=dist from cam
	float angle=10;float dist=0.6;
	
	//make rotation mvt
	HomogeneousMatrix posec_c2;
	posec_c2.TranslateZ(-dist);
	posec_c2.RotateY(-angle);
	posec_c2.TranslateZ(dist);
	
	//posec_c2.TranslateX(-dist*0.1);
	//posec_c2.TranslateZ(-dist*0.01);
	//posec_c2.TranslateZ(-dist*0.1);
	//posec_c2.RotateY(1.);
	
	std::cout<<"translationGT= "<<sqrt(posec_c2.get_translation().squaredNorm())<<std::endl;
	
	//create Points that will do measures
	int nb_point_all=4000;	
	//keep projections of points saved to be able to have same noise
	Vector2f projectionPoint[nb_point_all][nb_kf];
	Vector2f projectionPointInPair[nb_point_all][nb_kf];
	
	Vector3f RandomPoints[nb_point_all];
	for(int i=0;i<nb_point_all;i++)
	{
		//get coord
		Vector3f pt;
		float angle2=2*3.1415*((double)rand()/(double)RAND_MAX);
		pt[0]=cos(angle2)*(dist*2);
		pt[1]=0;
		pt[2]=sin(angle2)*(dist*2);
		
		for(int j=0;j<3;j++)pt[j]+=0.2*((double)rand()/(double)RAND_MAX-0.5);
		//std::cout<<((double)rand()/(double)RAND_MAX-0.5)<<std::endl;
		RandomPoints[i]=pt;
	}	
	
	
	HomogeneousMatrix poseTemp;
	//poseTemp.Init(0.01,0.01,0.01,0,0,0);
	poseTemp.TranslateZ(dist);
	poseTemp.RotateZ(1.57);
	cv::Mat img(im_h,im_w,CV_8UC1);
	img.setTo(0);
	
	//create all Keyframes
	for(int i=0;i<360/angle+10 && i<nb_kf;i++)
	{	
		Map_Estim.createNewKeyFrame(img,poseTemp);	
		PoseGT[i]=poseTemp;
		KeyFrame &kfi=*Map_Estim.getKF(i);
		
		//create best features
		//directly without creating bestMatches => not possible to do standard BA
		/*int i1p_c=0;
		for(int j=0;j<nb_point_all;j++)
		{
			//try project point and if good then create best feature local
			Vector3f ptCam=poseTemp*RandomPoints[j];
			if(ptCam[2]>0)
			{
				Vector2f proj=myCamera.Project(ptCam);
				Vector2f noise;
				noise[0]=gaussianNoise()*image_noise;
				noise[1]=gaussianNoise()*image_noise;
				proj+=noise;
				
				if(proj[0]>0 && proj[0]<im_w && proj[1]>0 && proj[1]<im_h)
				{
					uptoscaleFeature newFeat;
					newFeat.posRef=myCamera.ToMeters(proj);//position of feature in meters in KF
					newFeat.depthInRef=ptCam[2];//estimated depth
					newFeat.recAngle=2.;//reconstruction angle	
					newFeat.scoreFundamentalOrigin=10.;
					newFeat.i1p=i1p_c;//id of corresponding corner
					i1p_c++;
		
					newFeat.nb_outlier=0;
		
	#ifdef SAVE_POINT_COLOR
					newFeat.col[0]=0;
					newFeat.col[1]=0;
					newFeat.col[2]=0;
	#endif
					kfi.addLocalBestFeature(newFeat);
				}
			}
		}*/
		
		//create best features using local good stereo
		HomogeneousMatrix poseMiniStereo;
		HomogeneousMatrix relMiniStereo(-0.05,0,0,0,0,0);
		poseMiniStereo=relMiniStereo*poseTemp;
		
		kfi.setBestRelPose(relMiniStereo);
		std::vector<p_match> bestMatches;
		
		
		int i1p_c=0;
		for(int j=0;j<nb_point_all;j++)
		{
			//try project point and if good then create best feature local
			Vector3f ptCam=poseTemp*RandomPoints[j];
			Vector3f ptCam2=poseMiniStereo*RandomPoints[j];
			//if(ptCam[2]>0)
			if(ptCam[2]>0 && ptCam2[2])
			{
				//std::cout<<j<<std::endl;
				Vector2f proj=myCamera.Project(ptCam);
				Vector2f noise;
				noise[0]=gaussianNoise()*image_noise;
				noise[1]=gaussianNoise()*image_noise;
				proj+=noise;
				
				Vector2f proj2=myCamera.Project(ptCam2);
				Vector2f noise2;
				noise2[0]=gaussianNoise()*image_noise;
				noise2[1]=gaussianNoise()*image_noise;
				proj2+=noise2;
				
				projectionPoint[j][i]=proj;
				projectionPointInPair[j][i]=proj2;
				
				if(proj[0]>0 && proj[0]<im_w && proj[1]>0 && proj[1]<im_h)
				if(proj2[0]>0 && proj2[0]<im_w && proj2[1]>0 && proj2[1]<im_h)
				{
					//std::cout<<j<<std::endl;
					uptoscaleFeature newFeat;
					newFeat.posRef=myCamera.ToMeters(proj);//position of feature in meters in KF
					newFeat.depthInRef=ptCam[2];//estimated depth
					newFeat.recAngle=2.;//reconstruction angle	
					newFeat.scoreFundamentalOrigin=10.;
					newFeat.i1p=i1p_c;//id of corresponding corner
		
					newFeat.nb_outlier=0;
		
	#ifdef SAVE_POINT_COLOR
					//newFeat.grayVal=0;
					for(int c=0;c<3;c++)
					newFeat.col[c]=0;
	#endif
					kfi.addLocalBestFeature(newFeat);
					
					//add match in best matches
					p_match newMatch(proj[0],proj[1],i1p_c,proj2[0],proj2[1],i1p_c);
					bestMatches.push_back(newMatch);
					i1p_c++;
				}
			}
		}
		kfi.setBestLocalMatches(bestMatches);
		
		//best stereo has been created using, want it to be set using algebraic solution
		HomogeneousMatrix relPoseNoscale;
		bool isMotionGood=CheckNewStereo(bestMatches,&myCamera,relPoseNoscale);
		if(isMotionGood)
		{
			kfi.clearLocalBestFeatures();
			
			//rescale to gt norm
			float normNoScale=sqrt(relPoseNoscale.get_translation().squaredNorm());
			float normGT=sqrt(relMiniStereo.get_translation().squaredNorm());
			relPoseNoscale.set_translation((normGT/normNoScale)*relPoseNoscale.get_translation());
			kfi.setBestRelPose(relPoseNoscale);
			
			
			for(int j=0;j<bestMatches.size();j++)
			{
				float depthInRef;
				float recAngle;
				
				Vector2f mes1=myCamera.ToMeters(Vector2f(bestMatches[j].u1p,bestMatches[j].v1p));//in ref
				Vector2f mes2=myCamera.ToMeters(Vector2f(bestMatches[j].u1c,bestMatches[j].v1c));//in current
				reconstructionFromRays(mes1,mes2,relPoseNoscale,depthInRef,recAngle,false);
				
				uptoscaleFeature newFeature;
				newFeature.posRef=myCamera.ToMeters(Vector2f(bestMatches[j].u1p,bestMatches[j].v1p));
				newFeature.depthInRef=depthInRef;
				newFeature.recAngle=recAngle;
				newFeature.i1p=bestMatches[j].i1p;
			
#ifdef SAVE_POINT_COLOR
				for(int c=0;c<3;c++)
					newFeature.col[c]=0;
#endif
			
				kfi.addLocalBestFeature(newFeature);
			}
			kfi.doMiniBA(&myCamera);
			
		}
		else
			coutRed<<"Matches no good for mini stereo init"<<endlRed;
		
		
		
		
		poseTemp=posec_c2*poseTemp;
	}
		
	//create neigboring
	std::cout<<"create neigboring"<<std::endl;
	/*for(int i=0;i<nb_kf-1;i++)
	{
		KeyFrame &kf1=*Map_Estim.getKF(i);	  
		KeyFrame &kf2=*Map_Estim.getKF(i+1);*/	 
	for(int c=0;c<2;c++)
	for(int i=0;i<nb_kf-1-c;i++)
	{
		KeyFrame &kf1=*Map_Estim.getKF(i);	  
		KeyFrame &kf2=*Map_Estim.getKF(i+c+1);	 
		std::cout<<"create edge between kf"<<i<<" and "<<i+1<<std::endl;
		
		HomogeneousMatrix relPose=kf1.getPose()*kf2.getPose().inverse();
		//kf2 is kfc and kf1 is kfn
		
		
		NeigbourKFNew newNeigbor1;
		newNeigbor1.neighboring_kf=kf2.getId();
		newNeigbor1.relative_poses=relPose;
		newNeigbor1.Homography=Matrix3f::Identity();//no tracking here so do not care

		MatrixXf InformationMatrixPose(6,6);InformationMatrixPose.setZero();
		
		//create matches
		std::vector<p_match> matches_2_to_1;
		int i1p_1=0;
		int i1p_2=0;
		for(int j=0;j<nb_point_all;j++)
		{
			//try project point and if good then create best feature local
			bool point_visible_in_1=false;
			Vector3f ptCam1=kf1.getPose()*RandomPoints[j];
			Vector3f ptCam1m=kf1.getBestRelPose()*kf1.getPose()*RandomPoints[j];
			Vector2f proj1;
			Vector2f proj1m;
			if(ptCam1[2]>0 && ptCam1m[2]>0)
			{
				/*proj1=myCamera.Project(ptCam1);
				Vector2f noise;
				noise[0]=gaussianNoise()*image_noise;
				noise[1]=gaussianNoise()*image_noise;
				proj1+=noise;
				
				proj1m=myCamera.Project(ptCam1m);
				Vector2f noise2;
				noise2[0]=gaussianNoise()*image_noise;
				noise2[1]=gaussianNoise()*image_noise;
				proj1m+=noise;*/
				proj1=projectionPoint[j][i];
				proj1m=projectionPointInPair[j][i];

				
				if(proj1[0]>0 && proj1[0]<im_w && proj1[1]>0 && proj1[1]<im_h)
				if(proj1m[0]>0 && proj1m[0]<im_w && proj1m[1]>0 && proj1m[1]<im_h)
					point_visible_in_1=true;

			}
			
			bool point_visible_in_2=false;
			Vector3f ptCam2=kf2.getPose()*RandomPoints[j];
			Vector3f ptCam2m=kf2.getBestRelPose()*kf2.getPose()*RandomPoints[j];
			Vector2f proj2;
			Vector2f proj2m;
			if(ptCam2[2]>0 && ptCam2m[2]>0)
			{

				proj2=projectionPoint[j][i+c+1];
				proj2m=projectionPointInPair[j][i+c+1];

				if(proj2[0]>0 && proj2[0]<im_w && proj2[1]>0 && proj2[1]<im_h)
				if(proj2m[0]>0 && proj2m[0]<im_w && proj2m[1]>0 && proj2m[1]<im_h)
					point_visible_in_2=true;
			}
			
			if(point_visible_in_1 && point_visible_in_2)
			{
				p_match mMatch(proj2[0],proj2[1],i1p_2,proj1[0],proj1[1],i1p_1);
				matches_2_to_1.push_back(mMatch);
				
				//compute information matrix
				int idFeat_c=kf1.indexCandidateFeatureFromVisoId(mMatch.i1c);
				int idFeat_n=kf2.indexCandidateFeatureFromVisoId(mMatch.i1p);
				if(idFeat_c!=-1 && idFeat_n!=-1)
				{
					
					uptoscaleFeature *feat_c=kf1.getPtLocalBestFeatures(idFeat_c);
					uptoscaleFeature *feat_n=kf2.getPtLocalBestFeatures(idFeat_n);

					Vector3f scaleError=feat_c->getLocalCoordinates()-(relPose*feat_n->getLocalCoordinates());
					//deriv with respect to relPose:
					MatrixXf jac_dp(3,6);
					jac_dp.block(0,0,3,3)=-Matrix3f::Identity();
					jac_dp.block(0,3,3,3)=GetSkew(relPose*(feat_n->getLocalCoordinates()));	

					InformationMatrixPose+=jac_dp.transpose()*jac_dp;		
				}
			}
			
			if(point_visible_in_1)i1p_1++;
			if(point_visible_in_2)i1p_2++;
			
		}
		std::cout<<"nb matches created "<<matches_2_to_1.size()<<std::endl;
		
		/*std::cout<<"InformationMatrixPoseC = "<<std::endl;
		std::cout<<InformationMatrixPoseC<<std::endl<<std::endl;
		std::cout<<"InformationMatrixPoseN = "<<std::endl;
		std::cout<<InformationMatrixPoseN<<std::endl<<std::endl;
		
		Eigen::FullPivLU<MatrixXf> lu(InformationMatrixPoseC);
		std::cout<<"rank InformationMatrixPoseC = "<<lu.rank()<<std::endl;
		MatrixXf tri_u = lu.matrixLU().triangularView<Upper>();
		std::cout<<"tri_uC = "<<std::endl;
		std::cout<<tri_u<<std::endl;
		
		Eigen::FullPivLU<MatrixXf> lun(InformationMatrixPoseC);
		std::cout<<"rank InformationMatrixPoseN = "<<lun.rank()<<std::endl;
		MatrixXf tri_un = lun.matrixLU().triangularView<Upper>();
		std::cout<<"tri_un = "<<std::endl;
		std::cout<<tri_un<<std::endl;*/
		
		//newNeigbor1.InformationMatrixPose=InformationMatrixPose;//Information matrix of current kf of newNeigbor1 = kf1 => correspond to InformationMatrixPoseN
		
		newNeigbor1.matches=matches_2_to_1;
		newNeigbor1.edgeScore=10;//will be updated if 3D points are generated by this link
		kf1.addNeighbour(newNeigbor1);
		
		NeigbourKFNew newNeigbor2;
		newNeigbor2.neighboring_kf=kf1.getId();
		newNeigbor2.relative_poses=kf2.getPose()*kf1.getPose().inverse();
		//newNeigbor2.InformationMatrixPose=InformationMatrixPose;
		
		newNeigbor2.Homography=Matrix3f::Identity();//no tracking here so do not care
		for(int m=0;m<matches_2_to_1.size();m++)
			newNeigbor2.matches.push_back(matches_2_to_1[m].reverse());
		

		newNeigbor2.edgeScore=10;//will be updated if 3D points are generated by this link
		kf2.addNeighbour(newNeigbor2);
		
		//need to use getRelativePoseAndScale after getNeigbor as the function uses the edges between the keyframes
		NeigbourKFNew &newNeigbor1ref=*kf1.getPtNeigbour(kf1.getNbNeigbours()-1);
		NeigbourKFNew &newNeigbor2ref=*kf2.getPtNeigbour(kf2.getNbNeigbours()-1);
		Map_Estim.getRelativePoseAndScale(kf1.getId(),kf2.getId(),newNeigbor1ref.relative_scale,newNeigbor1ref.relative_poses,newNeigbor1ref.Informationscale,newNeigbor1ref.InformationMatrixPose);
		std::cout<<"scale constraint opt neigb kf"<<kf1.getId()<<"= "<<newNeigbor1ref.relative_scale<<std::endl;
		Map_Estim.getRelativePoseAndScale(kf2.getId(),kf1.getId(),newNeigbor2ref.relative_scale,newNeigbor2ref.relative_poses,newNeigbor2ref.Informationscale,newNeigbor2ref.InformationMatrixPose);
		std::cout<<"scale constraint opt neigb kf"<<kf2.getId()<<"= "<<newNeigbor2ref.relative_scale<<std::endl;
		
	}
	
	//create map points
	std::cout<<"create LUT"<<std::endl;
	int LUTtoPtId[nb_point_all];
	int LUTtoKForigin[nb_point_all];
	for(int i=0;i<nb_point_all;i++)
	{
		LUTtoPtId[i]=-1;
		LUTtoKForigin[i]=-1;
	}
	//int PtPerKf[nb_kf];
	//for(int i=0;i<nb_kf;i++)PtPerKf[i]=0;
	
	for(int j=0;j<nb_point_all;j++)
	{
		//check number of projections
		int nb_proj=0;
		for(int k=0;k<nb_kf;k++)
		{
			KeyFrame &kfc=*Map_Estim.getKF(k);
			Vector3f ptCam=kfc.getPose()*RandomPoints[j];
			Vector3f ptCam2=kfc.getBestRelPose()* kfc.getPose()*RandomPoints[j];
			if(ptCam[2]>0 && ptCam2[2]>0)
			{
				//Vector2f proj=myCamera.Project(ptCam);
				//Vector2f proj2=myCamera.Project(ptCam2);
				Vector2f proj=projectionPoint[j][k];
				Vector2f proj2=projectionPointInPair[j][k];
				if(proj[0]>0 && proj[0]<im_w && proj[1]>0 && proj[1]<im_h)
				if(proj2[0]>0 && proj2[0]<im_w && proj2[1]>0 && proj2[1]<im_h)
				{
					nb_proj++;
					if(nb_proj==2)
					{
						int currentNbPtInKf=kfc.getNbMapPoint();
						LUTtoPtId[j]=currentNbPtInKf;
						LUTtoKForigin[j]=k;
						
						MapPoint _newPoint;
						_newPoint.setId(currentNbPtInKf);
						_newPoint.updatePosition(RandomPoints[j]);
		#ifdef SAVE_POINT_COLOR			
						//_newPoint.setGrayVal(0);
						unsigned char col[3];
						for(int c=0;c<3;c++)col[c]=0;
						_newPoint.setCol(col);

		#endif
						_newPoint.setWeight(10);
				
						//add point to KF
						//add it to the neigbor
						kfc.addMapPoint(_newPoint);
					}
				}
			}
		}
	}
	
	
	std::cout<<"use LUT to create links points to features"<<std::endl;
	//create views and matches //ie links between points and features
	for(int i=0;i<nb_kf;i++)
	{
		KeyFrame &kfc=*Map_Estim.getKF(i);	  
		
		int i1p_c=0;
		for(int j=0;j<nb_point_all;j++)
		{
			//try project point and if good then create best feature local
			Vector3f ptCam=kfc.getPose()*RandomPoints[j];
			Vector3f ptCam2=kfc.getBestRelPose()* kfc.getPose()*RandomPoints[j];
			if(ptCam[2]>0 && ptCam2[2]>0)
			{
				//Vector2f proj=myCamera.Project(ptCam);
				//Vector2f proj2=myCamera.Project(ptCam2);
				Vector2f proj=projectionPoint[j][i];
				Vector2f proj2=projectionPointInPair[j][i];
				if(proj[0]>0 && proj[0]<im_w && proj[1]>0 && proj[1]<im_h)
				if(proj2[0]>0 && proj2[0]<im_w && proj2[1]>0 && proj2[1]<im_h)
				{
					//point is view from 
					if(LUTtoPtId[j]!=-1)
					{
						KeyFrame *KfOrigin=Map_Estim.getKF(LUTtoKForigin[j]);
						int idPtInKf=LUTtoPtId[j];
						MapPoint &point=*KfOrigin->getPtMapPoint(idPtInKf);
						
						//add match
						
						int id_feat_c=kfc.indexCandidateFeatureFromVisoId(i1p_c);//id_feat_c should be equal to i1p_c
						uptoscaleFeature &featc= *kfc.getPtLocalBestFeatures(id_feat_c);
						featc.matched=true;
						featc.idPoint=idPtInKf;
						featc.ptKForigin=KfOrigin;
						
						//add view
						point.addView(i,id_feat_c);
					}
					
					
					i1p_c++;
				}
				
			}
		}
	}
	
	//for the moment point have GT position; now use local features in views to update estimated position
	for(int i=0;i<nb_kf;i++)
	{
		KeyFrame &KFc=*Map_Estim.getKF(i);	
		for(int p=0;p<KFc.getNbMapPoint();p++)
		{
			
			//if(p<10)std::cout<<"\t\tpoint "<<p<<std::endl;
			//now we have the point and each of its views
			//each view has to be used to compute the sum to minimise
			MapPoint *point=KFc.getPtMapPoint(p);
			
			float newWeight=0;
			Vector3f newPosition=Vector3f(0,0,0);
			if(point->isUsed())//if point is good
			for(int v=0;v<point->nbViews();v++)
			{
				//if(p<10)std::cout<<"\t\t\tview = "<<v<<std::endl;
				int kfView=point->getView(v);
				//if(p<10)std::cout<<"\t\t\tkfview = "<<kfView<<std::endl;
				int i1pView=point->getI1p(v);
				KeyFrame &KFv=*Map_Estim.getKF(kfView);
				//get corresponding local feature:
				int id_feat_v=KFv.indexCandidateFeatureFromVisoId(i1pView);
				//if(p<10)std::cout<<"\t\t\tid_feat_v = "<<id_feat_v<<std::endl;
				if(id_feat_v!=-1)
				{
					uptoscaleFeature &feat_v=*KFv.getPtLocalBestFeatures(id_feat_v);
					
					newWeight+=feat_v.scoreFundamentalOrigin;							
					newPosition+=feat_v.scoreFundamentalOrigin*(KFv.getPose().inverse()*feat_v.getLocalCoordinates());
				}
			}
			point->updatePosition(newPosition/newWeight);
			point->setWeight(newWeight);
		}				
				
	}
	
}

//draw map viewed from kf_disp if disp_best_local ==0
//else draw map from best ministereo pair of kf_disp
void addDrawFunction(void) 
{	
	KeyFrame &kfd=*Map_Estim.getKF(kf_disp);
  
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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
#include <iostream>
#include <fstream>
using namespace std;

void optimiseAndOutputForPlotting()
{
	int nbItertor=100;
	ofstream myfile;
	//set init position
	for(int i=0;i<2;i++)
		disturb();
	
	//save init position
	/*HomogeneousMatrix pose0[nb_kf];
	for(int k=0;k<nb_kf;k++)
	      pose0[k]=Map_Estim.getKF(k)->getPose();*/
	Map_Estim.saveToFile("mapTest.map");
      
	//set init pose and do optim
	//for(int k=0;k<nb_kf;k++)
	//      Map_Estim.getKF(k)->setPose(pose0[k]);
	Map_Estim.loadFromFile("mapTest.map");

	
	//should do BA at end since modify inner structure of Keyframes => save us some resetting to be done
	myfile.open ("ReprojEvolEss.txt");
	float minEss;
	for(int iter=0;iter<nbItertor;iter++)
	{
		float reprojErr=Map_Estim.getReprojectionError();
		myfile << iter<<"\t"<<reprojErr <<"\n";
		doMapOptimEssential();
		
		if(iter==0)
		      minEss=reprojErr;
		else
		      if(reprojErr<minEss)minEss=reprojErr;
	}
	myfile.close();
	
	//set init pose and do optim
	Map_Estim.loadFromFile("mapTest.map");
	
	//should do BA at end since modify inner structure of Keyframes => save us some resetting to be done
	myfile.open ("ReprojEvolVar.txt");
	float minVar;
	for(int iter=0;iter<nbItertor;iter++)
	{
		float reprojErr=Map_Estim.getReprojectionError();
		myfile << iter<<"\t"<<reprojErr <<"\n";
		doClassicMinFeatureVariance();	
		if(iter==0)
		      minVar=reprojErr;
		else
		      if(reprojErr<minVar)minVar=reprojErr;
	}
	myfile.close();
	
	//set init pose and do optim
	Map_Estim.loadFromFile("mapTest.map");
	
	//should do BA at end since modify inner structure of Keyframes => save us some resetting to be done
	myfile.open ("ReprojEvolVar2.txt");
	float minVar2;
	for(int iter=0;iter<nbItertor;iter++)
	{
		float reprojErr=Map_Estim.getReprojectionError();
		myfile << iter<<"\t"<<reprojErr <<"\n";
		doClassicMinFeatureVariance2();	
		if(iter==0)
		      minVar2=reprojErr;
		else
		      if(reprojErr<minVar2)minVar2=reprojErr;
	}
	myfile.close();
	
	//set init pose and do optim
	Map_Estim.loadFromFile("mapTest.map");
	
	//should do BA at end since modify inner structure of Keyframes => save us some resetting to be done
	myfile.open ("ReprojEvolBA.txt");
	float minBA;
	for(int iter=0;iter<nbItertor;iter++)
	{
		float reprojErr=Map_Estim.getReprojectionError();
		myfile << iter<<"\t"<<reprojErr <<"\n";
		doClassicBA();
		if(iter==0)
		      minBA=reprojErr;
		else
		      if(reprojErr<minBA)minBA=reprojErr;
	}
	myfile.close();
	
	std::cout<<"min errors:"<<std::endl;
	std::cout<<"BA:"<<minBA<<std::endl;
	std::cout<<"Var:"<<minVar<<std::endl;
	std::cout<<"Var2:"<<minVar2<<std::endl;
	std::cout<<"Ess:"<<minEss<<std::endl;
}
