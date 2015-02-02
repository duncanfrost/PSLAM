#include "MapTracker.h"
#include "../VisuOdo/motionEstim.h"

MapTracker::MapTracker(Camera *_cam,obMap *_map)
{
	Init(_cam,_map);
}

void MapTracker::Init(Camera *_cam,obMap *_map)
{
	myCamera=_cam;
	myMap=_map;	
	idrelKF=-1;
	MotionPriorHomography=Matrix3f::Identity();	
}


void MapTracker::TrackFrame(cv::Mat &_current_img,cv::Mat *_current_img_col,bool useCol)
{
    //std::cout<<"#############################################"<<std::endl;


    //coutBlue << "Number of keyframes: " << myMap->getNbKeyFrames() << endlBlue;
	//start lony knowing previous position with respect to closest KF in map
	//ie the one with an overlap of imahe space > KFoverlap and closest in term of translation
	
	//prepare current image pymid
	cv::Mat img_c[NB_LEVELS];
	_current_img.copyTo(img_c[0]);
	for(int i=1;i<NB_LEVELS;i++)
		cv::pyrDown( img_c[i-1], img_c[i], cv::Size( img_c[i-1].cols/2, img_c[i-1].rows/2 ) );	


	//=> try to do matching with each of the KF
	//use the one with the most matches to do coarse pose estimation
	idrelKF=0;
	int nb_match_best=-1;
	amoTimer timerMatching;
	timerMatching.start();

	Matrix3f newMotionPriorHomography;
#ifndef USE_OMP
    //std::cout<<"MapTracker : update each close KF with new frame"<<std::endl;
	for(int n=0;n<id_closestKF.size();n++)
	{
        //std::cout<<"\tKF["<<id_closestKF[n]<<"] useNewFrame "<<std::endl;
		
		Matrix3f currentHomog=myMap->getKF(id_closestKF[n])->getHomography();
		Matrix3f currentHomogPlusPrior=MotionPriorHomography*currentHomog;
		myMap->getKF(id_closestKF[n])->setHomography(currentHomogPlusPrior);
		int nb_matchn=myMap->getKF(id_closestKF[n])->useNewFrame(img_c,myCamera);
		
		if(nb_match_best<nb_matchn)
		{
			nb_match_best=nb_matchn;
			idrelKF=id_closestKF[n];
			//if best match then use it to compute new displacement prior
			Eigen::FullPivLU<MatrixXf> lu(currentHomog);		
			newMotionPriorHomography=myMap->getKF(id_closestKF[n])->getHomography()*lu.inverse();
		}

	}
#else
	omp_lock_t lock;
	omp_init_lock(&lock);
	
	std::cout<<"MapTracker : update each close KF with new frame"<<std::endl;
	#pragma omp parallel for
	for(int n=0;n<id_closestKF.size();n++)
	{
		omp_set_lock(&lock);										
		std::cout<<"\tKF["<<id_closestKF[n]<<"] useNewFrame "<<std::endl;
		omp_unset_lock(&lock);
		
		//bit to compute motion prior
		Matrix3f currentHomog=myMap->getKF(id_closestKF[n])->getHomography();
		Matrix3f currentHomogPlusPrior=MotionPriorHomography*currentHomog;		
		//myMap->getKF(id_closestKF[n])->setHomography(currentHomogPlusPrior);//to turn displacement prior on/off just comment or uncomment this line
		
		int nb_matchn=myMap->getKF(id_closestKF[n])->useNewFrame(img_c,myCamera);
		
		//myMap->getKF(id_closestKF[n])->computeOverlap();
		//std::cout<<"overlap after track = "<<myMap->getKF(id_closestKF[n])->getOverlapWithLastFrame()<<std::endl;
		
		omp_set_lock(&lock);										
		if(nb_match_best<nb_matchn)
		{
			nb_match_best=nb_matchn;
			idrelKF=id_closestKF[n];
			//if best match then use it to compute new displacement prior
			Eigen::FullPivLU<MatrixXf> lu(currentHomog);		
			newMotionPriorHomography=myMap->getKF(id_closestKF[n])->getHomography()*lu.inverse();
			
		}
		omp_unset_lock(&lock);
	}
	omp_destroy_lock(&lock);
#endif
	
	if(nb_match_best!=-1)
		MotionPriorHomography=newMotionPriorHomography;
	
	timerMatching.stop("Matching with KFs");
    //std::cout<<"nb_match_best = "<<nb_match_best<<std::endl;
	
	
	if(nb_match_best>MIN_MATCHES_EDGE)	
	{
        //std::cout<<"MapTracker : checkNeigborsForPoints"<<std::endl;
		amoTimer timerLinkage;
		timerLinkage.start();
		//introduction of currentFrame in one neigbor KF could have given it better local stereo
		//if that s the case then, if KF linked to map then new points and views might be
		//introduced. => check among neigbors if we need to check for new points
		for(int n=0;n<id_closestKF.size();n++)
			if(myMap->getKF(id_closestKF[n])->haveLocalBestFeaturesChanged())
				myMap->checkNeigborsForPoints(id_closestKF[n]);
			
		timerLinkage.stop("Linkage");
		
		//use best neigboring KF to refine pose estimation
		//for that use all map points and local features linked to matches
        //std::cout<<"pose refine before Optim map"<<std::endl;
		KeyFrame &closestKF=*myMap->getKF(idrelKF);
		//relPose=closestKF.computeRelativeCurrentPoseWithAllMatches(myCamera);
		
		//at that stage have new current frame and might have new matches and points in map
		//=> do very local Bundle adjustment around current position
        //std::cout<<"optimiseInnerWindow"<<std::endl;
		
		//get KF in inner window
		int depth_inner=2;
		std::vector<int> optimKF;
		optimKF=myMap->getConnectedKeyframes(id_closestKF,depth_inner);
		
		int nb_iter=5;bool robust =true;
		//can use minimisation variance local features
		MapOptimiser mBA(myMap);
		//mBA.optimiseInnerWindowRobust(optimKF,nb_iter);//use only matches=>could be easy to marginalised but misses some information
		mBA.optimiseInnerWindowRobust2(optimKF,nb_iter);//use mapPoints that gather all views, uses all info but hard to marginalise
		
		//can use minimisation of epipolar geometry cost and min variance features
		//MapOptimiserEssential mBA(myMap);
		//mBA.optimiseInnerWindowRobust(optimKF,nb_iter);
		
        //std::cout<<"pose refine after Optim map"<<std::endl;
		//now use best 3D estimation of current KF and matches with current image
		//to get relative pose;
		
		//can get relative pose from what has been computed already in keyframe
		//Indeed: in each keyframe have fundamental matrix decomposed and ministereo reconstructed
		//this is rescaled to fit features of Keyframe
		//however from one fundamental matrix decomposition to an other noise has big effect=> give very coarse result
		//relPose=closestKF.getRelativePose();

		//prefer estimating pose by minimixing reprojection error between measures and keyframe best matched
		//if last frame is best fundamental then computeRelativeCurrentPose should not change estimation of pose...
		
		//can use only local features => pb can have jumps ..
		//relPose=closestKF.computeRelativeCurrentPoseWithLocalFeatures(myCamera);
		
		
		//pose with points from BA
		//use map points: pb need to have the local map optimised first
		//relPose=closestKF.computeRelativeCurrentPoseWithMatchedFeatures(myCamera);
		relPose=closestKF.computeRelativeCurrentPoseWithAllMatches(myCamera);
		
		
		
	
		amoTimer timerKeyFrameProcess;
		timerKeyFrameProcess.start();
		
		//check if need new Keyframe:
		//do if max overlap with closest KFs is less than limitoverlap
		float max_overlap=0;
		for(int n=0;n<id_closestKF.size();n++)
		{
			KeyFrame &fk=*myMap->getKF(id_closestKF[n]);
			//get overlap
			float current_overlap=fk.getOverlapWithLastFrame();
			std::cout<<"current_overlap KF["<<id_closestKF[n]<<"]: "<<current_overlap<<std::endl;
			if(current_overlap>max_overlap)
				max_overlap=current_overlap;
		}
        //std::cout<<"max overlap: "<<max_overlap<<std::endl;
		
        if((max_overlap<KFoverlap) && (myMap->getNbKeyFrames() < 1))
		{
			coutRed<<"small overlap: "<<max_overlap<<", create new KF"<<endlRed;
			if(!useCol)
				myMap->createNewKeyFrame(img_c[0],relPose,idrelKF,id_closestKF);
			else
				myMap->createNewKeyFrame(*_current_img_col,relPose,idrelKF,id_closestKF);
			
			//std::vector<int> optimKF2;
			//optimKF2=myMap->getConnectedKeyframes(id_closestKF);			
			//mBA.optimiseInnerWindowRobust(optimKF2,nb_iter);
		}
		else
		{
			//check min reconstruction angle with all keyframe having overlap < min_overlap
			//if angle big enough then create keyframe
			float min_recangle=3.14/2;//set as max value rec angle could have
			for(int n=0;n<id_closestKF.size();n++)
			{
				KeyFrame &fk=*myMap->getKF(id_closestKF[n]);
				//get overlap
				float current_overlap=fk.getOverlapWithLastFrame();
				if(current_overlap>KFoverlap)
				{
					float recangle=fk.getAverageReconstructionAngle();
					//std::cout<<"recangle : "<<recangle<<std::endl;
					if(recangle<min_recangle)
						min_recangle=recangle;
				}
			}
			//do if reconstruction angle is large enough
            if(min_recangle>KFreconstructionAngle && myMap->getNbKeyFrames() < 1)
			{
				
				/*cv::imwrite("ref1.png",_current_img);//for offline
				std::cout<<"Homography = "<<std::endl;
				std::cout<<myMap->getKF(0)->getHomography()<<std::endl;
				
				exit(1);*/
				//coutRed<<"WARNING rec angle has been change for debug"<<endlRed;
				coutGreen<<"good reconstruction angle: "<<min_recangle<<", create new KF"<<endlGreen;
				if(!useCol)
					myMap->createNewKeyFrame(img_c[0],relPose,idrelKF,id_closestKF);
				else
					myMap->createNewKeyFrame(*_current_img_col,relPose,idrelKF,id_closestKF);
				
				//std::vector<int> optimKF2;
				//optimKF2=myMap->getConnectedKeyframes(id_closestKF);			
				//mBA.optimiseInnerWindowRobust(optimKF2,nb_iter);

			}
		}
		
		
		//update close KFs		
		//remove ones with too small overlap
		for(int n=0;n<id_closestKF.size();n++)
		{
			KeyFrame &fk=*myMap->getKF(id_closestKF[n]);
			myMap->removeUnusedPoints(id_closestKF[n]);
			
			//get overlap
			float current_overlap=fk.getOverlapWithLastFrame();
			if(current_overlap<KFoverlapMin)
			{
				//myMap->removeUnusedPoints(id_closestKF[n]);
				id_closestKF.erase(id_closestKF.begin()+n);
				n--;
			}
			
			//update their last relative pose
			HomogeneousMatrix currentRelPose=getPose()*fk.getPose().inverse();
			fk.setRelativePosePrevious(currentRelPose);
			
		}
		

        //std::cout<<"id_closestKF = "<<std::endl;
//		for(int i=0;i<id_closestKF.size();i++)std::cout<<id_closestKF[i]<<"  ";
//		std::cout<<std::endl;
		while(id_closestKF.size()>MAX_KF_ACTIVE)
		{
			std::vector<int>::iterator it =  std::min_element(id_closestKF.begin(), id_closestKF.end());
			id_closestKF.erase(it);
			
		}
        //std::cout<<"id_closestKF after = "<<std::endl;
        //for(int i=0;i<id_closestKF.size();i++)std::cout<<id_closestKF[i]<<"  ";
        //std::cout<<std::endl;

		timerKeyFrameProcess.stop("KeyFrame Process");
		
		//for(int k=0;k<myMap->getNbKeyFrames();k++)
		//	myMap->checkMapPointIds(k);
	}
	else
	{
		coutRed<<"lost => create new submap"<<endlRed;
		
		//else we were lost => need new Keyframe
		int idNewKF;
		if(!useCol)
			idNewKF=myMap->createUnconnectedNewKF(img_c[0]);
		else
			idNewKF=myMap->createUnconnectedNewKF(*_current_img_col);
		
		//remove previous connection to map
		id_closestKF.clear();
		//rotation_vs_KF.clear();
		
		//create new one
		id_closestKF.push_back(idNewKF);
		//rotation_vs_KF.push_back(Vector3f(0,0,0));//not rotation between current and new KF => 0
		
		//init relPose as identity
		relPose=HomogeneousMatrix();
		MotionPriorHomography=Matrix3f::Identity();
		
		//cv::imwrite("ref0.png",_current_img);
	}
	
	//for(int n=0;n<myMap->getNbKeyFrames();n++)
	//	coutRed<<"KF["<<n<<"] nb mpa point = "<<myMap->getKF(n)->getNbMapPoint()<<endlRed;
	
	

}

