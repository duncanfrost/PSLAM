#include "KeyFrame.h"

//#define USE_MINIBA

KeyFrame::KeyFrame()
{
	InitMemory();
}

void KeyFrame::InitMemory()
{
	matchesCurrent.reserve(700);
  
	best_matches.reserve(700);
	best_features.reserve(1000);
	
	neigbours.reserve(10);
	mMapPoints.reserve(600);
	localBestFeaturesHaveChanged=false;
	overlapWithLast=1;
}

KeyFrame::KeyFrame(int _id,cv::Mat &_img,HomogeneousMatrix22 _pose)
{
	InitMemory();
	Init(_id,_img,_pose);
}

void KeyFrame::Init(int _id,cv::Mat &_img,HomogeneousMatrix22 _pose)
{
	w_To_cam=_pose;
	id=_id;
	cv::Mat imgBW;//in case have color img as param
	cv::Mat *imgBWpt=NULL;
	
	if(_img.type()!=CV_8UC1)
	{
		_img.copyTo(img_col);
		cv::cvtColor(_img,imgBW,CV_RGB2GRAY);
		imgBWpt=&imgBW;
	}
	else
	{
		cv::cvtColor(_img,img_col,CV_GRAY2RGB);
		imgBWpt=&_img;
	}

	makeKF_Image(*imgBWpt);
	extractORBFeatures();
	initMatcher();

	 	
	//set default depth range
	depthRange[0]=0.1;
	depthRange[1]=1.;
	
	//kf not matched with anything yet:
	best_score_fundamental=0;
	//scale=1;
	Homography=Matrix3f::Identity();
	//MotionPriorHomography=Matrix3f::Identity();
	average_recAngle=0;
	last_fundamental_score=0;
	img_best_pair=imgBWpt->clone();
}

void KeyFrame::makeKF_Image(cv::Mat &_img)
{
	_img.copyTo(img_p[0]);
	for(int i=1;i<NB_LEVELS;i++)
		cv::pyrDown( img_p[i-1], img_p[i], cv::Size( img_p[i-1].cols/2, img_p[i-1].rows/2 ) );	
}



void KeyFrame::initMatcher()
{
	matcher.InitWithRef(img_p);
}
KeyFrame::~KeyFrame()
{

}


bool CheckNewStereo(std::vector<p_match> &matches,Camera *myCamera,HomogeneousMatrix22 &KFtoCurrent)
{
	if(matches.size()>10)
	{
		//get up to scale 3D displacement with all matches (slow)
		MotionEstimation mMotionEstimator(myCamera->getFocalPix()[0],myCamera->getFocalPix()[1],myCamera->getCenterImage()[0],myCamera->getCenterImage()[1]);
		bool isEstimationValid=mMotionEstimator.estimateMotion(matches);
		//bool isEstimationValid=mMotionEstimator.estimateMotionAndRemoveOutliers(matches);
		
		if(isEstimationValid)
		{
			amoTimer timerCheckFundamental;
			timerCheckFundamental.start();
			
			//std::cout<<"KFtoCurrent computation"<<std::endl;
			//check how is 3d reconstruction from this mini stereo
			Matrix3f Rotation;Vector3f translation;
			mMotionEstimator.getEstimatedTransfo(&Rotation(0,0),&translation[0]);
			KFtoCurrent.Init(translation,Rotation);//default translation normalized to 1
			return true;
		}  
	}
	return false;
}

std::vector<uptoscaleFeature> KeyFrame::getGoodFeaturesFromRayIntersection(std::vector<p_match> &matches,Camera *myCamera,HomogeneousMatrix22 &KFtoCurrent)
{
	std::vector<uptoscaleFeature> GoodLocalFeatures;
	for(int i=0;i<matches.size();i++)
		if(matches[i].i1p!=-1)//can only happen when init local stereo with best neigbor: only then i1p can be = -1
		  //since matches are defined with inverse matches with best neigbor. If ==-1 then all the info of the point provide
		//is unsed in estimation of fundamental. Later will not be able to match it with other feature so useless to create local feature with it
	{
		float depthInRef;
		float recAngle;
		
		Vector2f mes1=myCamera->ToMeters(Vector2f(matches[i].u1p,matches[i].v1p));//in ref
		Vector2f mes2=myCamera->ToMeters(Vector2f(matches[i].u1c,matches[i].v1c));//in current
		
		reconstructionFromRays(mes1,mes2,KFtoCurrent,depthInRef,recAngle,true);
		//reconstructionFromRays(mes1,mes2,KFtoCurrent,depthInRef,recAngle,false);//can turn of intersection verification but that adds noise in mapping...
		
		if(depthInRef>0 && !isinf(depthInRef))//reconstruction valid
		{
			uptoscaleFeature newFeature;
			newFeature.posRef=myCamera->ToMeters(Vector2f(matches[i].u1p,matches[i].v1p));
			newFeature.depthInRef=depthInRef;
			newFeature.recAngle=recAngle;
			newFeature.i1p=matches[i].i1p;
			
#ifdef SAVE_POINT_COLOR
			//newFeature.grayVal=img_p[0].at<unsigned char>(matches[i].v1p,matches[i].u1p);
			newFeature.col[0]=img_col.at<cv::Vec3b>(matches[i].v1p,matches[i].u1p)[2];
			newFeature.col[1]=img_col.at<cv::Vec3b>(matches[i].v1p,matches[i].u1p)[1];
			newFeature.col[2]=img_col.at<cv::Vec3b>(matches[i].v1p,matches[i].u1p)[0];
#endif
			
			GoodLocalFeatures.push_back(newFeature);
		}
		else
		{
			//if want to use miniBA then need to remove outlier matches for now
#ifdef USE_MINIBA
			matches.erase(matches.begin()+i);
			i--;
#endif
		}
	}
	return GoodLocalFeatures;
}
std::vector<uptoscaleFeature> KeyFrame::filterWithDepthConsistancy(std::vector<uptoscaleFeature> &GoodLocalFeatures)
{
	std::vector<uptoscaleFeature> FilteredLocalFeatures;
	if(GoodLocalFeatures.size()>10)
	{
		//get depthRange
		float minDepth=GoodLocalFeatures[0].depthInRef;
		float maxDepth=GoodLocalFeatures[0].depthInRef;
		for(int i=0;i<GoodLocalFeatures.size();i++)
		{
			if(GoodLocalFeatures[i].depthInRef<minDepth)minDepth=GoodLocalFeatures[i].depthInRef;
			if(GoodLocalFeatures[i].depthInRef>maxDepth)maxDepth=GoodLocalFeatures[i].depthInRef;
		}
		float depthRange=maxDepth-minDepth;
		
		//check that first neigbor have similar depth, if not don t keep it
		if(depthRange==0)
			FilteredLocalFeatures=GoodLocalFeatures;
		else
		for(int i=0;i<GoodLocalFeatures.size();i++)
		{
			//find closest feature
			int id_closest=0;
			if(i==0)id_closest=1;
			
			float dist_closest=(GoodLocalFeatures[i].posRef-GoodLocalFeatures[id_closest].posRef).transpose()*(GoodLocalFeatures[i].posRef-GoodLocalFeatures[id_closest].posRef);
			
			for(int j=1;j<GoodLocalFeatures.size();j++)
				if(i!=j)
				{
					float dist_j=(GoodLocalFeatures[i].posRef-GoodLocalFeatures[j].posRef).transpose()*(GoodLocalFeatures[i].posRef-GoodLocalFeatures[j].posRef);
					if(dist_closest>dist_j)
					{
						id_closest=j;
						dist_closest=dist_j;
					}
				}
				
			//check that depth difference is inferior to percentage of depthRange
			float diffDepth=(GoodLocalFeatures[i].depthInRef-GoodLocalFeatures[id_closest].depthInRef);
			diffDepth=(diffDepth<0)?-diffDepth:diffDepth;
			
			if(diffDepth/depthRange<0.1)
			{
				FilteredLocalFeatures.push_back(GoodLocalFeatures[i]);
			}
		}
	}
	return FilteredLocalFeatures;
}
float KeyFrame::getFundamentalMatrixScore(std::vector<uptoscaleFeature> &GoodFeatures,float &_avr_angle)
{
	//accumulate reconstruction angle of final valid feature to check fundamental matrix score
	//=> best score if we have lots of features with large reconstruction angle
	//advantage is that it is scale invariant
	float accumulationRecAngle=0;
	for(int i=0;i<GoodFeatures.size();i++)accumulationRecAngle+=GoodFeatures[i].recAngle;
	
	//std::cout<<"Fundamental matrix score   : "<<accumulationRecAngle<<std::endl;
	float score_fund=0;
	if(GoodFeatures.size()!=0)
		score_fund=accumulationRecAngle;
	if(GoodFeatures.size()!=0)
		_avr_angle=accumulationRecAngle/GoodFeatures.size();
	else
		_avr_angle=0;
	
	return score_fund;
}

struct MeasureTempPoint
{
	Vector2f mMeasure;//position in meter in image
	Vector3f mPosition;//position of 3D point linked to measure
	float weight;
	unsigned char fromPoint;
	short KFsrc;
};

//predicate to use find if and features
struct CompareId
{
  CompareId(int Id) : id_(Id) {}
  bool operator()(const uptoscaleFeature feat) const {
    return id_ == feat.i1p;
  }
  private:
    int id_;
};

int KeyFrame::indexCandidateFeatureFromVisoId(int idFeaturep)
{
	std::vector<uptoscaleFeature>::iterator it = std::find_if (best_features.begin(), best_features.end(), CompareId(idFeaturep));
	if(it != best_features.end())
		return it-best_features.begin();
	else 
		return -1;
	
}

HomogeneousMatrix22 KeyFrame::computeRelativeCurrentPoseWithLocalFeatures(Camera *_myCamera)
{
	//get list of MeasureTempPointfrom matches and best features
	std::vector<MeasureTempPoint> measurePoints;
	for(int m=0;m<matchesCurrent.size();m++)
	{
		int idFeaturep=matchesCurrent[m].i1p;
		//check if we got a 3D position estimation for the corresponding point
		std::vector<uptoscaleFeature>::iterator it = std::find_if (best_features.begin(), best_features.end(), CompareId(idFeaturep));
		if(it != best_features.end())
		{
			//match has a 3D position linked to its feature in reference
			MeasureTempPoint newMeasure;
			newMeasure.mMeasure=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1c,matchesCurrent[m].v1c));
			
			//measure in ref
			Vector2f mMeasureRef=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1p,matchesCurrent[m].v1p));
			newMeasure.mPosition=toHomogeneous(mMeasureRef)* it->depthInRef;
			measurePoints.push_back(newMeasure);
		}
		
	}
	
	if(measurePoints.size()>10)
	{
		//estimate relative pose by minimizing reprojection error:
		int nb_iter=20;
		
		for(int iter=0;iter<nb_iter;iter++)
		{

			VectorXf Jte(6);Jte.setZero();
			MatrixXf H(6,6);H.setZero();
			
			
			for(int i=0;i<measurePoints.size();i++)
			{
				//matrix to be filled			
				MeasureTempPoint measure=measurePoints[i];
				Vector3f mapPointsCam=relPose*measure.mPosition;
				
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=measure.mMeasure;//desired projection= measurement
				Vector2f x_c=_myCamera->ProjectZ1(mapPointsCam);//current projection
				
				Vector2f error=x_d-x_c;
				
				float norm_reproj_error=error.transpose()*error;

				//get jacobien  of error with respect to variation of camera pose
				MatrixXf de_dp=-_myCamera->ProjectZ1_Jac_Dp(mapPointsCam);	
				
				Jte+=de_dp.transpose()*error;
				
				H+=de_dp.transpose()*de_dp;
				
			}
			Eigen::FullPivLU<MatrixXf> lu(H);
			//todo try p+=Dp
			VectorXf Dp(6);
			Dp=-1.*(lu.inverse()*Jte);
            relPose=HomogeneousMatrix22(Dp)* relPose;
		}
	}
	return relPose;
}
HomogeneousMatrix22 KeyFrame::computeRelativeCurrentPoseWithMatchedFeatures(Camera *_myCamera)
{
	//get list of MeasureTempPointfrom matches and best features
	std::vector<MeasureTempPoint> measurePoints;
	for(int m=0;m<matchesCurrent.size();m++)
	{
		int idFeaturep=matchesCurrent[m].i1p;
		//check if we got a 3D position estimation for the corresponding point
		std::vector<uptoscaleFeature>::iterator it = std::find_if (best_features.begin(), best_features.end(), CompareId(idFeaturep));
		if(it != best_features.end())
		{
			if(it->matched)
			{
				//match has a 3D position linked to its feature in reference
				MeasureTempPoint newMeasure;
				newMeasure.mMeasure=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1c,matchesCurrent[m].v1c));
				
				//measure in ref
				Vector2f mMeasureRef=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1p,matchesCurrent[m].v1p));
				newMeasure.mPosition=w_To_cam*it->ptKForigin->getPtMapPoint(it->idPoint)->getPosition();
				measurePoints.push_back(newMeasure);
			}
		}
		
	}
	
	if(measurePoints.size()>10)
	{
		//estimate relative pose by minimizing reprojection error:
		int nb_iter=20;
		
		for(int iter=0;iter<nb_iter;iter++)
		{

			VectorXf Jte(6);Jte.setZero();
			MatrixXf H(6,6);H.setZero();
			
			
			for(int i=0;i<measurePoints.size();i++)
			{
				//matrix to be filled			
				MeasureTempPoint measure=measurePoints[i];
				Vector3f mapPointsCam=relPose*measure.mPosition;
				
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=measure.mMeasure;//desired projection= measurement
				Vector2f x_c=_myCamera->ProjectZ1(mapPointsCam);//current projection
				
				Vector2f error=x_d-x_c;
				
				float norm_reproj_error=error.transpose()*error;

				//get jacobien  of error with respect to variation of camera pose
				MatrixXf de_dp=-_myCamera->ProjectZ1_Jac_Dp(mapPointsCam);	
				
				Jte+=de_dp.transpose()*error;
				
				H+=de_dp.transpose()*de_dp;
				
			}
			Eigen::FullPivLU<MatrixXf> lu(H);
			//todo try p+=Dp
			VectorXf Dp(6);
			Dp=-1.*(lu.inverse()*Jte);
            relPose=HomogeneousMatrix22(Dp)* relPose;
		}
	}
	return relPose;

}
HomogeneousMatrix22 KeyFrame::computeRelativeCurrentPoseWithAllMatches(Camera *_myCamera)
{
	//get list of MeasureTempPointfrom matches and best features
	std::vector<MeasureTempPoint> measurePoints;
	std::cout<<"collect measures "<<matchesCurrent.size()<<" matches "<<std::endl;
	int matchToMapPoint=0;
	int matchToFeat=0;
	for(int m=0;m<matchesCurrent.size();m++)
	{
		int idFeaturep=matchesCurrent[m].i1p;
		//check if we got a 3D position estimation for the corresponding feature
		std::vector<uptoscaleFeature>::iterator it = std::find_if (best_features.begin(), best_features.end(), CompareId(idFeaturep));
		if(it != best_features.end())
		{
		  
			if(it->matched && it->ptKForigin->getPtMapPoint(it->idPoint)->isUsed())//logically there should be no feature matched to unsued point
			{
				//match has a map point linked to its feature in reference
				MeasureTempPoint newMeasure;
				newMeasure.mMeasure=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1c,matchesCurrent[m].v1c));
				
				//measure in ref
				Vector2f mMeasureRef=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1p,matchesCurrent[m].v1p));
				MapPoint &point=*it->ptKForigin->getPtMapPoint(it->idPoint);
				newMeasure.mPosition=w_To_cam* point.getPosition();
				
				newMeasure.weight=point.getWeight();
				//newMeasure.weight=1.;
				
				newMeasure.fromPoint=1;//for testing
				newMeasure.KFsrc=it->ptKForigin->id;
				measurePoints.push_back(newMeasure);
				
				matchToMapPoint++;
			}
			else
			{
				//match has a loca feature linked to its feature in reference
				MeasureTempPoint newMeasure;
				newMeasure.mMeasure=_myCamera->ToMeters(Vector2f(matchesCurrent[m].u1c,matchesCurrent[m].v1c));
				
				//measure in ref
				Vector2f mMeasureRef=it->posRef;
				newMeasure.mPosition=toHomogeneous(mMeasureRef)* it->depthInRef;
				
				Vector2f measPix=_myCamera->ToPixels(mMeasureRef);
				//if(it->i1p==41)
				//	std::cout<<"create measure "<<it->i1p<<": "<<matchesCurrent[m].u1p<<"  "<<matchesCurrent[m].v1p<<" <=>" <<matchesCurrent[m].u1c<<"  "<<matchesCurrent[m].v1c<<" <=> "<<measPix[0]<<" "<<measPix[1]<<std::endl;
				
				newMeasure.weight=it->scoreFundamentalOrigin;
				//newMeasure.weight=1.;
				
				newMeasure.fromPoint=0;//for testing
				measurePoints.push_back(newMeasure); 
				
				matchToFeat++;
			}
		}
		
	}
	//std::cout<<"matchToMapPoint = "<<matchToMapPoint<<std::endl;
	//std::cout<<"matchToFeat = "<<matchToFeat<<std::endl;
	
	if(measurePoints.size()>10)
	{
		//estimate relative pose by minimizing reprojection error:
		int nb_iter=30;
		float mLMLambda = 0.0001;//before LevMarConstant
		float mdLambdaFactor = 2.0;
		
		for(int iter=0;iter<nb_iter;iter++)
		{
			//get tukey factor for robust estimation
			std::vector<float> vdErrorSquared;
			for(int m=0;m<measurePoints.size();m++)
			{
				MeasureTempPoint measure=measurePoints[m];
				Vector3f mapPointsCam=relPose*measure.mPosition;
				if(mapPointsCam[2]>0)
				{
					//compute error with observed points (in meter in z=1 plane)
					Vector2f x_d=measure.mMeasure;//desired projection= measurement
					Vector2f x_c=_myCamera->ProjectZ1(mapPointsCam);//current projection
					
					Vector2f errorPix=_myCamera->m2PixProjJac()*(x_d-x_c);
					//if(measure.fromPoint==0)std::cout<<"\tLocal: "<<_myCamera->ToPixels(x_d).transpose()<<"  \tcurrent proj: "<<_myCamera->ToPixels(x_c).transpose()<<std::endl;
					//if(measure.fromPoint==1)std::cout<<"\tPoint: "<<_myCamera->ToPixels(x_d).transpose()<<"  \tcurrent proj: "<<_myCamera->ToPixels(x_c).transpose()<<std::endl;
					vdErrorSquared.push_back(errorPix.transpose()*errorPix);
				}
			}

			//std::cout<<"\tvdErrorSquared.size() "<<vdErrorSquared.size()<<std::endl;
			
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
			
			int nbf=0;
			float residuef=0;
			int nbp=0;
			float residuep=0;
			
			float residue=0;
			float weightTotal=0;
			float max_residue=0;
			int nb_meas_used=0;
			for(int i=0;i<measurePoints.size();i++)
			{
				//matrix to be filled			
				MeasureTempPoint measure=measurePoints[i];
				Vector3f mapPointsCam=relPose*measure.mPosition;
				if(mapPointsCam[2]>0)
				{
					//compute error with observed points (in meter in z=1 plane)
					Vector2f x_d=measure.mMeasure;//desired projection= measurement
					Vector2f x_c=_myCamera->ProjectZ1(mapPointsCam);//current projection
					
					Vector2f error=x_d-x_c;
					float norm_reproj_error=(_myCamera->m2PixProjJac()*error).squaredNorm();
					float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);					
					if(TukeyCoef>0)
					{
						float norm_reproj_error=error.transpose()*error;						

						//get jacobien  of error with respect to variation of camera pose
						MatrixXf de_dp=-_myCamera->ProjectZ1_Jac_Dp(mapPointsCam);	
						
						Jte+=measure.weight* de_dp.transpose()*error;
						
						H+=measure.weight*de_dp.transpose()*de_dp;
						
						float err_pix=sqrt((_myCamera->m2PixProjJac()*error).squaredNorm());
						if(err_pix>max_residue)max_residue=err_pix;
						residue+=TukeyCoef*err_pix;
						weightTotal+=TukeyCoef;
						if(measure.fromPoint==0){residuef+=err_pix;nbf++;}
						if(measure.fromPoint==1){residuep+=err_pix;nbp++;}
						nb_meas_used++;
					}
				}
				
			}
			//std::cout<<"\tweightTotal "<<weightTotal<<std::endl;
			
			//std::cout<<"\tmax_residue = "<<max_residue<<std::endl;
			/*if(iter == 0)
			{
				std::cout<<"\tresidue = "<<residue/measurePoints.size()<<std::endl;
				std::cout<<"\tresiduef = "<<residuef/nbf<<std::endl;
				std::cout<<"\tresiduep = "<<residuep/nbp<<std::endl;				  
			}*/
			if(nb_meas_used<4)break;
			
			for(int i=0;i<6;i++)
				H(i,i)=(1.+mLMLambda)*H(i,i);
			
			Eigen::FullPivLU<MatrixXf> lu(H);
			//todo try p+=Dp
			VectorXf Dp(6);
			Dp=-1.*(lu.inverse()*Jte);
			if(!isnan(Dp[0]))
			{
			
                HomogeneousMatrix22 relPoseNew=HomogeneousMatrix22(Dp)* relPose;
				
				float residueAfter=0;
				float weightTotalAfter=0;
				for(int i=0;i<measurePoints.size();i++)
				{
					//matrix to be filled			
					MeasureTempPoint measure=measurePoints[i];
					Vector3f mapPointsCam=relPoseNew*measure.mPosition;
					if(mapPointsCam[2]>0)
					{
						//compute error with observed points (in meter in z=1 plane)
						Vector2f x_d=measure.mMeasure;//desired projection= measurement
						Vector2f x_c=_myCamera->ProjectZ1(mapPointsCam);//current projection
						
						Vector2f error=x_d-x_c;
						float norm_reproj_error=(_myCamera->m2PixProjJac()*error).squaredNorm();
						float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);					
						if(TukeyCoef>0)
						{
							float err_pix=sqrt((_myCamera->m2PixProjJac()*error).squaredNorm());
							weightTotalAfter+=TukeyCoef;
							residueAfter+=TukeyCoef*err_pix;
						}
					}
					
				}
				//std::cout<<"\tweightTotalAfter "<<weightTotalAfter<<std::endl;
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
	return relPose;

}
HomogeneousMatrix22 KeyFrame::RescalePose(HomogeneousMatrix22 &relPoseNoscale,std::vector<uptoscaleFeature> &NewFeatures,std::vector<uptoscaleFeature> &FilteredNewFeatures,std::vector<uptoscaleFeature> &best_features)
{
	//rescale relPoseNoscale and features FilteredNewFeatures to match best_features
	//std::cout<<"FilteredNewFeatures.size() = "<<FilteredNewFeatures.size()<<std::endl;
	std::vector<float> DesiredSurCurrent;
	for(int i=0;i<FilteredNewFeatures.size();i++)
	{
		//check if same feature has its depth estimated in best_features
		int idFeature=FilteredNewFeatures[i].i1p;
		//check if we got a 3D position estimation for the corresponding point
		std::vector<uptoscaleFeature>::iterator it = std::find_if (best_features.begin(), best_features.end(), CompareId(idFeature));
		if(it!=best_features.end())	
		{
			//std::cout<<"Depth desired = "<<it->depthInRef<<std::endl;
			//std::cout<<"Current desired = "<<FilteredNewFeatures[i].depthInRef<<std::endl;
			float scale_t=it->depthInRef/FilteredNewFeatures[i].depthInRef;
			if(!isnan(scale_t))
				DesiredSurCurrent.push_back(scale_t);
		}
	}
		
	if(DesiredSurCurrent.size()>0)
	{
		//have set of current to desired depth to estimate scale
		//do median filter (could do ransac too but median filter seems to work well)
		std::sort (DesiredSurCurrent.begin(), DesiredSurCurrent.end());
		float best_scale_ransac=DesiredSurCurrent[DesiredSurCurrent.size()/2];
		
		//refine using robust approach	
		float CurrentToDesiredScale;
		if(DesiredSurCurrent.size()==1)
			CurrentToDesiredScale=best_scale_ransac;
		else
		{
			std::vector<float> vdErrorSquared;
			for(int i=0;i<DesiredSurCurrent.size();i++)
			{
				float error=best_scale_ransac-DesiredSurCurrent[i];
				vdErrorSquared.push_back(error*error);
			}
			float sigma_tukey=getSigmaSquared(vdErrorSquared);
			
			float totScale=0;
			float totWeight=0;
			for(int i=0;i<DesiredSurCurrent.size();i++)
			{
				float error=best_scale_ransac-DesiredSurCurrent[i];
				float rootErr=sqrt(error*error);
				float TukeyCoef=squareRootTukey(rootErr,sigma_tukey);
				if(TukeyCoef!=0)
				{
					totScale+=TukeyCoef*DesiredSurCurrent[i];
					totWeight+=TukeyCoef;
				}
			}
			CurrentToDesiredScale=totScale/totWeight;
		}
		
		std::cout<<"\tkf["<<id<<"] DesiredSurCurrent.size() = "<<DesiredSurCurrent.size()<<std::endl;	
		std::cout<<"\tkf["<<id<<"] CurrentToDesiredScale = "<<CurrentToDesiredScale<<std::endl;	
		
		if(CurrentToDesiredScale<=0)//failed
		  return relPoseNoscale;
		  
		//std::cout<<"best_scale_ransac = "<<best_scale_ransac<<std::endl;
		
		//rescale features and transformation
		for(int i=0;i<NewFeatures.size();i++)
			NewFeatures[i].depthInRef=CurrentToDesiredScale*NewFeatures[i].depthInRef;
		for(int i=0;i<FilteredNewFeatures.size();i++)
			FilteredNewFeatures[i].depthInRef=CurrentToDesiredScale*FilteredNewFeatures[i].depthInRef;
		
        HomogeneousMatrix22 KFtoCurrent=relPoseNoscale;
		KFtoCurrent.set_translation(CurrentToDesiredScale*KFtoCurrent.get_translation());
		return KFtoCurrent;
	}
	else 
	{
		//don t know scale, most likely because there has not been enough translation => set t = 0
        HomogeneousMatrix22 KFtoCurrent=relPoseNoscale;
		KFtoCurrent.set_translation(Vector3f(0,0,0));
		return KFtoCurrent;
	}
	
}

int KeyFrame::useNewFrame(cv::Mat *_img_c,Camera *_myCamera)
{
  
	average_recAngle=0;
	last_fundamental_score=0;
	
	amoTimer timerMatching;
	//timerMatching.start();
	matcher.match(_img_c);
	//timerMatching.stop("matching");
	Homography=matcher.getHomography();
	computeOverlap();//update overlap with ref
	matchesCurrent=matcher.getCurrentMatches();
	int nb_matches=matchesCurrent.size();
	
	std::cout<<"\tKF["<<id<<"] nb matchesCurrent = "<<nb_matches<<std::endl;
	std::cout<<"\tKF["<<id<<"] overlap = "<<overlapWithLast<<std::endl;
	
	//get 3D if we can
	//get up to scale Homogeneous matrix from KF to current

    HomogeneousMatrix22 relPoseNoscale;
	bool isMotionGood=CheckNewStereo(matchesCurrent,_myCamera,relPoseNoscale);
	
	if(isMotionGood)
	{
	  
		//if manage to decompose fundamental then can try reconstructing features in 3D
		std::vector<uptoscaleFeature> NewFeatures=getGoodFeaturesFromRayIntersection(matchesCurrent,_myCamera,relPoseNoscale);
		std::cout<<"\tKF["<<id<<"]: nb NewFeatures = "<<NewFeatures.size()<<std::endl;		
		
		//do a bit of non linear optimisation with matches and relPoseNoscale
#ifdef USE_MINIBA
		doMiniBA(matchesCurrent,NewFeatures,_myCamera,relPoseNoscale);
#endif
		
		//can filter this 3D points to make it more robust: consider that depth should be consistent over the image
		//std::vector<uptoscaleFeature> FilteredNewFeatures=NewFeatures;
		std::vector<uptoscaleFeature> FilteredNewFeatures=filterWithDepthConsistancy(NewFeatures);
		//std::cout<<"\tKF["<<id<<"]: nb FilteredNewFeatures = "<<FilteredNewFeatures.size()<<std::endl;		
		
		//all the features are still with unknown scale
		//=> if any features in FilteredNewFeatures corresponds to the best ones stored in KF
		//then we can do rescale
		//At this stage if there are already map point linked to the KF then best_features should be aligned to them
		relPose=RescalePose(relPoseNoscale,NewFeatures,FilteredNewFeatures,best_features);
		//std::cout<<"\tKF["<<id<<"]: relPose = "<<relPose<<std::endl;
			  
		
		//get score from fundamental matrix
		last_fundamental_score=getFundamentalMatrixScore(FilteredNewFeatures,average_recAngle);
		std::cout<<"\tlast_fundamental_score = "<<last_fundamental_score<<std::endl;
		
		//check if that s the best stereo we had here
		if(best_score_fundamental<last_fundamental_score)
		{
			//SaveNewMatchesAndPoseForMiniBA(matchesCurrent,relPose,_myCamera);
			//check if we are linked to a map 
			if(getNbNeigbours()==0)
			{
				//must be first keyframe after getting lost
				//set average depth to default value
				float default_av_depth=1;
				float av_depth=0;
				
				for(int i=0;i<FilteredNewFeatures.size();i++)
					av_depth+=FilteredNewFeatures[i].depthInRef;
				av_depth=av_depth/FilteredNewFeatures.size();
				
				for(int i=0;i<NewFeatures.size();i++)
					NewFeatures[i].depthInRef=NewFeatures[i].depthInRef*default_av_depth/av_depth;
				for(int i=0;i<FilteredNewFeatures.size();i++)
					FilteredNewFeatures[i].depthInRef=FilteredNewFeatures[i].depthInRef*default_av_depth/av_depth;
				relPose.set_translation(relPose.get_translation()*default_av_depth/av_depth);
				
				
				img_best_pair=_img_c[0].clone();
				bestRelPose=relPose;
				best_score_fundamental=last_fundamental_score;
				best_matches=matchesCurrent;
				
				for(int i=0;i<NewFeatures.size();i++)
					NewFeatures[i].scoreFundamentalOrigin=best_score_fundamental;
				best_features=NewFeatures;

			}
			else
			{
				//will have to keep our existing good features since they are used in BA
				//=> add the new GoodFeatures that should already be rescaled to existing ones
				//some features are going to correspond to same corner=> in that case feature need to be replaced but need to keep matching information
				
				for(int i=0;i<NewFeatures.size();i++)
				{
					
					int idFeat=indexCandidateFeatureFromVisoId(NewFeatures[i].i1p);
					if(idFeat!=-1)//remove feature with same i1p
					{
						//copy matching information
						NewFeatures[i].matched=best_features[idFeat].matched;
						if(NewFeatures[i].matched)
						{
							NewFeatures[i].ptKForigin=best_features[idFeat].ptKForigin;
							NewFeatures[i].idPoint=best_features[idFeat].idPoint;
						}
						best_features.erase(best_features.begin()+idFeat);
					}

					NewFeatures[i].scoreFundamentalOrigin=best_score_fundamental;
					best_features.push_back(NewFeatures[i]);
				}
				
				img_best_pair=_img_c[0].clone();
				bestRelPose=relPose;
				best_matches=matchesCurrent;
				best_score_fundamental=last_fundamental_score;
				
				//will have to update mapPoints from new links that can be now matched with neigbors
				localBestFeaturesHaveChanged=true;
			}
		  
			//std::cout<<"KeyFrame::useNewFrame best_score_fundamental"<<std::endl;
			std::cout<<"\tKF["<<id<<"]: best_score_fundamental, best = "<<best_score_fundamental<<" current = "<<last_fundamental_score<<std::endl;
			std::cout<<"\tKF["<<id<<"]: best_features size = "<<NewFeatures.size()<<std::endl;
			//if it is then keep matches and score
			//but first need to check if can have scale information from neigbors
			//if no neigbors then no information on scale at all=> scale 1 is good by default
			
			//TODO but map points should be create while adding new KF
			//at that moment can check neigboring and tracks to create points

		}	
		else
		{
		  std::cout<<"\tKF["<<id<<"]: not best fundam, best = "<<best_score_fundamental<<" current = "<<last_fundamental_score<<std::endl;
		  //DoMiniBA(_myCamera);
		}

		//test if could do very small baselin BA

		
	}
	else
	{
	  	std::cout<<"\tKF["<<id<<"]: motion is not good, todo: could at least estimate rotation ?"<<std::endl;
	}

	return nb_matches;
}
/*void KeyFrame::displayMatchc(int i)
{
	for(int i=0;i<matchesCurrent.size();i++)
	{
		if(matchesCurrent[i].i1c==41)std::cout<<"featurec "<<matchesCurrent[i].i1c<<" : "<<matchesCurrent[i].u1c<<" "<<matchesCurrent[i].v1c<<std::endl;  
	}
}*/

void KeyFrame::useNeigborForInitLocalStereo(cv::Mat &imgBest,std::vector<p_match> &matches, HomogeneousMatrix22 &_relPoseNeigb,Camera *_myCamera)
{
	img_best_pair=imgBest.clone();
  
	average_recAngle=0;
	last_fundamental_score=0;
	
	matchesCurrent.clear();//KF should be all new so should be alread cleared.
	for(int i=0;i<matches.size();i++)
	{		
		//if(matches[i].i1c==41)std::cout<<"featurec "<<matches[i].i1c<<" : "<<matches[i].u1c<<" "<<matches[i].v1c<<std::endl;
		matchesCurrent.push_back(matches[i].reverse());
	}
	
	//_relPoseNeigb does not correspond to relative pose of keyframe with last current image !
	//it is the relative pose of this keyframe wrt its best neigbor=> not = to relPose
	bestRelPose=_relPoseNeigb.inverse();

	//get the 3D featrues from ray intersection should be already same scale as best neigbor
	best_features=getGoodFeaturesFromRayIntersection(matchesCurrent,_myCamera,bestRelPose);
	std::cout<<"\tnb best_features = "<<best_features.size()<<std::endl;
	//for(int i=0;i<10;i++)
	//	std::cout<<_myCamera->ToPixels(best_features[i].posRef).transpose()<<std::endl;

		
	//can filter this 3D points to make it more robust: consider that depth should be consistent over the image
	std::vector<uptoscaleFeature> FilteredNewFeatures=filterWithDepthConsistancy(best_features);
	//std::vector<uptoscaleFeature> FilteredNewFeatures=best_features;
		
	//get score from fundamental matrix
	best_score_fundamental=getFundamentalMatrixScore(FilteredNewFeatures,average_recAngle);
	for(int i=0;i<best_features.size();i++)
		best_features[i].scoreFundamentalOrigin=best_score_fundamental;
		
	std::cout<<"\tbest_score_fundamental, best = "<<best_score_fundamental<<std::endl;

	best_matches=matchesCurrent;
	matchesCurrent.clear();

}
void KeyFrame::extractORBFeatures()
{
	//pyramidal
	for(int l=0;l<NB_LEVELS;l++)
	{
		//ORB corners and descriptors
		if(l==1)//compute them only in first level of pyramid
		{
		std::vector<cv::KeyPoint> keypointsORB;
		cv::OrbFeatureDetector detector;
		//CV_WRAP explicit ORB(int nfeatures = 500, float scaleFactor = 1.2f, int nlevels = 8, int edgeThreshold = 31,
		//int firstLevel = 0, int WTA_K=2, int scoreType=ORB::HARRIS_SCORE, int patchSize=31 );	
		detector.detect( img_p[l], keypointsORB );
		for(int i = 0; i < keypointsORB.size(); i++)
			OrbCorners[l].push_back(Vector2f(keypointsORB[i].pt.x,keypointsORB[i].pt.y));
		
		//compute descriptors:
		//ORB		
		cv::OrbDescriptorExtractor extractor;
		extractor.compute( img_p[l], keypointsORB, OrbDescriptors[l] );		
		if(OrbDescriptors[l].type()!=CV_32F)
			OrbDescriptors[l].convertTo(OrbDescriptors[l], CV_32F);
		}

	}
}



int KeyFrame::isNeigbour(int idnewNeighbr)
{
	for(int i=0;i<neigbours.size();i++)
		if(neigbours[i].neighboring_kf==idnewNeighbr)
		{
			return i;
		}
	return 	-1;
}

int KeyFrame::getNumberOfFeatureMatched()
{
	int res=0;
	for(int i=0;i<best_features.size();i++)if(best_features[i].matched)res++;
	return res;
}
void KeyFrame::computeOverlap()
{
	//need to check warp and back warp: indeed for instance if zoom out, all ref is in current, but only part of current is in ref=>new info to get
	float overlap=getOverlapFromHomography(img_p[0].size(),Homography,4);
	Eigen::FullPivLU<MatrixXf> lu(Homography);		
	Matrix3f _Hinv=lu.inverse();	
	float overlapinv=getOverlapFromHomography(img_p[0].size(),_Hinv,4);
	overlapWithLast= (overlap<overlapinv)?overlap:overlapinv;
	
}

	//io functions
void KeyFrame::saveToStream(std::ofstream &fout, KeyFrame *kf0)
{
	fout.write((const char*)&id,sizeof(int));
	w_To_cam.saveToStream(fout);
	//std::cout<<"w_To_cam "<<w_To_cam<<std::endl;
	
	for(int i=0;i<3;i++)for(int j=0;j<3;j++)
	      fout.write((const char*)&Homography(i,j),sizeof(float));
	
	int _w=img_p[0].size().width;
	int _h=img_p[0].size().height;
	fout.write((const char*)&_w,sizeof(int));
	fout.write((const char*)&_h,sizeof(int));
	fout.write((const char*)(&img_p[0].at<unsigned char>(0,0)),img_p[0].size().width*img_p[0].size().height*sizeof(char));

	//depth range
	for(int i=0;i<2;i++)
	      fout.write((const char*)&depthRange[i],sizeof(float));
	
	fout.write((const char*)(&img_best_pair.at<unsigned char>(0,0)),img_best_pair.size().width*img_best_pair.size().height*sizeof(char));
	fout.write((const char*)&best_score_fundamental,sizeof(float));
	bestRelPose.saveToStream(fout);
	
	int nb_best_matches=best_matches.size();	
	fout.write((const char*)&nb_best_matches,sizeof(int));
	for(int i=0;i<nb_best_matches;i++)
		fout.write((const char*)&best_matches[i],sizeof(p_match));
		//best_matches[i].saveToStream(fout);

	int nb_features=best_features.size();	
	fout.write((const char*)&nb_features,sizeof(int));
	for(int i=0;i<nb_features;i++)
		best_features[i].saveToStream(fout,kf0);

	int nb_neigbours=neigbours.size();	
	fout.write((const char*)&nb_neigbours,sizeof(int));
	for(int i=0;i<nb_neigbours;i++)
		neigbours[i].saveToStream(fout);
	
	int nb_MapPoints=mMapPoints.size();	
	fout.write((const char*)&nb_MapPoints,sizeof(int));
	for(int i=0;i<nb_MapPoints;i++)
		mMapPoints[i].saveToStream(fout);
		
	
	//could save last info with current frame too...
	/*float average_recAngle;	
	float last_fundamental_score;
	
	//if had bestFeatures linked to map and have new ones=> have to check with neigbors if new connections are possible
	bool localBestFeaturesHaveChanged; 
	
	//lastely estimated relative pose between current frame and Keyframe
	HomogeneousMatrix relPose;*/
}
void KeyFrame::loadFromStream(std::ifstream &fout, KeyFrame *kf0)
{
	fout.read((char*)&id,sizeof(int));
	w_To_cam.loadFromStream(fout);

	for(int i=0;i<3;i++)for(int j=0;j<3;j++)
	      fout.read((char*)&Homography(i,j),sizeof(float));
	
	int width,height;	
	fout.read((char*)&width,sizeof(int));	
	fout.read((char*)&height,sizeof(int));
	cv::Mat img;	img.create(height,width, CV_8UC1);
	fout.read((char*)(&img.at<unsigned char>(0,0)),img.size().width*img.size().height*sizeof(char));
	makeKF_Image(img);
	cv::cvtColor(img,img_col,CV_GRAY2RGB);
	
	extractORBFeatures();
	initMatcher();	
	
	for(int i=0;i<2;i++)
	      fout.read((char*)&depthRange[i],sizeof(float));

	img_best_pair.create(height,width, CV_8UC1);
	fout.read((char*)(&img_best_pair.at<unsigned char>(0,0)),img_best_pair.size().width*img_best_pair.size().height*sizeof(char));
	fout.read((char*)&best_score_fundamental,sizeof(float));
	bestRelPose.loadFromStream(fout);
	
	int nb_best_matches;	
	fout.read((char*)&nb_best_matches,sizeof(int));
	best_matches.resize(nb_best_matches);
	for(int i=0;i<nb_best_matches;i++)
		fout.read((char*)&best_matches[i],sizeof(p_match));
		//best_matches[i].loadFromStream(fout);

	int nb_features;	
	fout.read((char*)&nb_features,sizeof(int));
	best_features.resize(nb_features);
	for(int i=0;i<nb_features;i++)
		best_features[i].loadFromStream(fout,kf0);
	
	int nb_neigbours;	
	fout.read((char*)&nb_neigbours,sizeof(int));
	neigbours.resize(nb_neigbours);
	for(int i=0;i<nb_neigbours;i++)
		neigbours[i].loadFromStream(fout);
	
	int nb_MapPoints;	
	fout.read((char*)&nb_MapPoints,sizeof(int));
	mMapPoints.resize(nb_MapPoints);
	for(int i=0;i<nb_MapPoints;i++)
	{
		mMapPoints[i].loadFromStream(fout);
		mMapPoints[i].setId(i);//should be unecessary
	}
}

void uptoscaleFeature::saveToStream(std::ofstream &fout, KeyFrame *kf0)
{
	for(int i=0;i<2;i++)
	      fout.write((const char*)&posRef[i],sizeof(float));		
	fout.write((const char*)&depthInRef,sizeof(float));		
	fout.write((const char*)&recAngle,sizeof(float));		
	fout.write((const char*)&scoreFundamentalOrigin,sizeof(float));		
	fout.write((const char*)&i1p,sizeof(int));		
	char ismatched=matched;
	fout.write((const char*)&ismatched,sizeof(char));
	int id_KFo=ptKForigin-kf0;
	fout.write((const char*)&id_KFo,sizeof(int));		
	fout.write((const char*)&idPoint,sizeof(int));		
#ifdef SAVE_POINT_COLOR
	//fout.write((const char*)&grayVal,sizeof(char));
	for(int i=0;i<3;i++)
		fout.write((const char*)&col[i],sizeof(char));
#endif
}
void uptoscaleFeature::loadFromStream(std::ifstream &fout, KeyFrame *kf0)
{
	for(int i=0;i<2;i++)
	      fout.read((char*)&posRef[i],sizeof(float));		
	fout.read((char*)&depthInRef,sizeof(float));		
	fout.read((char*)&recAngle,sizeof(float));		
	fout.read((char*)&scoreFundamentalOrigin,sizeof(float));		
	fout.read((char*)&i1p,sizeof(int));		
	char ismatched;
	fout.read((char*)&ismatched,sizeof(char));
	matched=(ismatched==1);
	int id_KFo;
	fout.read((char*)&id_KFo,sizeof(int));	
	ptKForigin=kf0+id_KFo;
	fout.read((char*)&idPoint,sizeof(int));			  
#ifdef SAVE_POINT_COLOR
	//fout.read((char*)&grayVal,sizeof(char));
	for(int i=0;i<3;i++)
		fout.read((char*)&col[i],sizeof(char));
#endif
  
}
//search for two unknown : alpha and beta so that alpha*v1-beta*v2=c1c2
//with alpha*v1 = c1p (p=intersection) and beta*v2=c2p
//=> c1p+pc2=c1c2
//express all the variables in frame1
int reconstructionFromRays(Vector2f mes1,Vector2f mes2,HomogeneousMatrix22 pose1_2,float &depth1, float &recAngle, bool checkIntersect)
{
    HomogeneousMatrix22 pose1;//identity
    HomogeneousMatrix22 pose2Inv=pose1_2.inverse();
	Vector3f centerCam2=pose2Inv.get_translation();

	Vector3f b=centerCam2;

	Vector3f pt_c1;pt_c1=toHomogeneous(mes1);
	Vector3f pt_c2;pt_c2=pose2Inv*toHomogeneous(mes2);

	Vector3f v1=pt_c1;
	Vector3f v2=pt_c2-centerCam2;

	MatrixXf A(3,2);
	for(int j=0;j<3;j++)A(j,0)=v1[j];
	for(int j=0;j<3;j++)A(j,1)=-v2[j];

	Matrix2f AtA;AtA=A.transpose()*A;
	
	if(AtA.determinant()==0)  //point at infinity
	{
		//std::cout<<"det = 0"<<std::endl;
		goto failIntersectionRays;

	}
	else
	{
		Eigen::FullPivLU<MatrixXf> lu(AtA);

		//x=(AtA)-1 At b
		Vector2f x=(lu.inverse()*A.transpose())*b;	
		
				//v1 and v2 have z=1 => x is really the depth and not distance to center
		float depth=x[0]*v1[2];
		float depth2=x[1]*v2[2];
		
		//check if depth > 0
		if(depth<=0 || depth2<=0)//point behind camera
		{
			//std::cout<<"depth < 0"<<std::endl;
			goto failIntersectionRays;
		}
		
		//check if two rays are actually intersecting
		//3d point from ray1:
		Vector3f X1=x[0]*v1;
		//3d point from ray2:
		Vector3f X2=centerCam2+x[1]*v2;
		//std::cout<<"diff = "<<(X1-X2).squaredNorm()/depth/sqrt(b.squaredNorm())<<std::endl;

		if(!checkIntersect || (X1-X2).squaredNorm()/depth/sqrt(b.squaredNorm())<1e-3 )
		{
			//success
			depth1=depth;
			//get reconstruction angle:
			v1.normalize();
			v2.normalize();
			recAngle=acos(v1.transpose()*v2);
			if(isnan(recAngle))recAngle=0;//can happen due to float things
			return 1;
		}
		else
			goto failIntersectionRays;
	}
	
	failIntersectionRays:
	depth1=-1;
	recAngle=0;
	return -1;
}

int reconstructionFromRaysInvDepth(Vector2f mes1,Vector2f mes2,HomogeneousMatrix22 pose1_2,float &invdepth, float &recAngle)
{
    HomogeneousMatrix22 pose1;//identity
    HomogeneousMatrix22 pose2Inv=pose1_2.inverse();
	Vector3f centerCam2=pose2Inv.get_translation();

	Vector3f b=centerCam2;

	Vector3f pt_c1;pt_c1=toHomogeneous(mes1);
	Vector3f pt_c2;pt_c2=pose2Inv*toHomogeneous(mes2);

	Vector3f v1=pt_c1;
	Vector3f v2=pt_c2-centerCam2;

	MatrixXf A(3,2);
	for(int j=0;j<3;j++)A(j,0)=v1[j];
	for(int j=0;j<3;j++)A(j,1)=-v2[j];

	Matrix2f AtA;AtA=A.transpose()*A;
	
	if(AtA.determinant()==0)  //point at infinity
	{
		invdepth=0;
		recAngle=0;
		return 1;

	}
	else
	{
		Eigen::FullPivLU<MatrixXf> lu(AtA);

		//x=(AtA)-1 At b
		Vector2f x=(lu.inverse()*A.transpose())*b;	
		
				//v1 and v2 have z=1 => x is really the depth and not distance to center
		float depth=x[0]*v1[2];
		float depth2=x[1]*v2[2];
		
		//check if depth > 0
		if(depth<=0 || depth2<=0)//point behind camera
		{
			//std::cout<<"depth < 0"<<std::endl;
			goto failIntersectionRays;
		}
		
		//check if two rays are actually intersecting
		//3d point from ray1:
		Vector3f X1=x[0]*v1;
		//3d point from ray2:
		Vector3f X2=centerCam2+x[1]*v2;
		//std::cout<<"diff = "<<(X1-X2).squaredNorm()/depth/sqrt(b.squaredNorm())<<std::endl;

		if( (X1-X2).squaredNorm()/depth/sqrt(b.squaredNorm())<1e-3 )
		{
			//success
			invdepth=1./depth;
			//get reconstruction angle:
			v1.normalize();
			v2.normalize();
			recAngle=acos(v1.transpose()*v2);
			if(isnan(recAngle))recAngle=0;//can happen due to float things
			return 1;
		}
		else
			goto failIntersectionRays;
	}
	
	failIntersectionRays:
	invdepth=-1;
	recAngle=0;
	return -1;
}

/*void KeyFrame::SaveNewMatchesAndPoseForMiniBA(std::vector<p_match> &matchesCurrent,HomogeneousMatrix &relPose,Camera *_myCamera)
{

	
	for(int m=0;m<nb_feat_ref;m++)
		listMatches[m][cpt_view_saved%max_view_in_mini_BA].matched=false;
	  
	for(int m=0;m<matchesCurrent.size();m++)
	{
		p_match &match=matchesCurrent[m];
		listMatches[match.i1p][cpt_view_saved%max_view_in_mini_BA].matched=true;
		listMatches[match.i1p][cpt_view_saved%max_view_in_mini_BA].pos_p=_myCamera->ToMeters(Vector2f(match.u1p,match.v1p));
		listMatches[match.i1p][cpt_view_saved%max_view_in_mini_BA].pos_c=_myCamera->ToMeters(Vector2f(match.u1c,match.v1c));
		
	}
	poseMiniBA[nb_view_saved%max_view_in_mini_BA]=relPose;
	
	cpt_view_saved++;
	if(nb_view_saved<max_view_in_mini_BA)
		nb_view_saved++;

}
struct miniBAjacz
{
      short index_pt_opt;//if max angle not big enough => pt not optim then index_pt_opt=-1
      short index_cam_opt;//if max angle not big enough => pt not optim then index_pt_opt=-1
      Vector2f proj_error;//deriv wrt point pos (only if index_pt_opt!=-1)
      Vector2f de_dz;//deriv wrt point pos (only if index_pt_opt!=-1)
      MatrixXf de_dp;//deriv wrt cam pos
};

void KeyFrame::DoMiniBA(Camera *myCamera)
{
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;
	if(nb_view_saved>0)
	for(int iter=0;iter<5;iter++)
	{
		//get invDepth
		for(int i=0;i<nb_feat_ref;i++)
		{
			//get corresponding 3d point
			float idepthp=invDepthFeat[i];
			if(idepthp==-1)//depth not initialised yet => use current match to initialise it
			{
				int idFeat=indexCandidateFeatureFromVisoId(i);
				if(idFeat!=-1)
				{
					idepthp=1./best_features[idFeat].depthInRef;
					invDepthFeat[i]=idepthp;
				}
			}
		}

	  
		//get Tukey factor
		std::cout<<"get Tukey factor iter "<<iter<<std::endl;
		float sigma_tukey[nb_view_saved];
		for(int k=0;k<nb_view_saved;k++)
		{		
			std::vector<float> vdErrorSquared;
			for(int i=0;i<nb_feat_ref;i++)
				if(listMatches[i][k].matched)
			{
				//get corresponding 3d point
				float idepthp=invDepthFeat[i];
				//if had good depth or if manage to initialise it, then can compute error
				
				
				if(idepthp!=-1)
				{
					if(idepthp<1e-2)idepthp=1e-2;
					//Vector3f mapPointsp=toHomogeneous(listMatches[i][k].pos_p)/idepthp;
					//Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(listMatches[i][k].pos_p)/idepthp);
					Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(MeanFeatPos[i])/idepthp);
					//std::cout<<"###########################################"<<std::endl;
					//std::cout<<"point in ref : "<<myCamera->ToPixels(listMatches[i][k].pos_p).transpose()<<std::endl;
					//std::cout<<"proj in ref : "<<myCamera->ToPixels(myCamera->ProjectZ1(mapPointsp)).transpose()<<std::endl;
					//std::cout<<"point in cur : "<<myCamera->ToPixels(listMatches[i][k].pos_c).transpose()<<std::endl;
					//std::cout<<"proj in cur : "<<myCamera->ToPixels(myCamera->ProjectZ1(mapPointsc)).transpose()<<std::endl;
					
					if(mapPointsc[2]>0)
					{
					    Vector2f x_c=myCamera->ProjectZ1(mapPointsc);//current projection
					    Vector2f error=myCamera->m2PixProjJac()*(listMatches[i][k].pos_c-x_c);//error in pixels
					    vdErrorSquared.push_back(error.transpose()*error);
					}
				}
			}
			std::cout<<"nb meas view ["<<k<<"] = "<<vdErrorSquared.size()<<std::endl;
			sigma_tukey[k]=getSigmaSquared(vdErrorSquared);
			if(sigma_tukey[k] < MINSIGMATUKEY)
				sigma_tukey[k] = MINSIGMATUKEY;	
		}
		
		
		//now collect all jacobians for camera and depth features to optimise (for that needs to be matched, have a positive depth which is not outlier)
		//get LUT for depths to optimise
		std::cout<<"get LUTs"<<std::endl;
		int LUT_to_opt[nb_feat_ref];
		int LUT_to_main[nb_feat_ref];
		
		int opt_cpt=0;
		for(int i=0;i<nb_feat_ref;i++)
		{
			LUT_to_opt[i]=-1;
			for(int k=0;k<nb_view_saved;k++)
			{
				if(listMatches[i][k].matched)
				{
					float idepthp=invDepthFeat[i];
					if(idepthp!=-1)
					{
						//Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(listMatches[i][k].pos_p)/idepthp);
						Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(MeanFeatPos[i])/idepthp);
						if(mapPointsc[2]>0)
						{
							Vector2f x_c=myCamera->ProjectZ1(mapPointsc);//current projection
							Vector2f error=myCamera->m2PixProjJac()*(listMatches[i][k].pos_c-x_c);//error in pixels
							float TukeyCoef=squareRootTukey(error.transpose()*error,sigma_tukey[k]);
							if(TukeyCoef>0)
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
		std::cout<<"get error and jacs"<<std::endl;
		float residue_robust=0;
		float valid_wmeas=0;
		std::vector<miniBAjacz> fJacobian;		

		for(int k=0;k<nb_view_saved;k++)
		{		
			for(int i=0;i<nb_feat_ref;i++)
				if(listMatches[i][k].matched)
			{
				//get corresponding 3d point
				float idepthp=invDepthFeat[i];
				if(idepthp!=-1)//depth not initialised yet => use current match to initialise it
				{
					//Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(listMatches[i][k].pos_p)/idepthp);
					Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(MeanFeatPos[i])/idepthp);
					if(mapPointsc[2]>0)
					{
						Vector2f x_c=myCamera->ProjectZ1(mapPointsc);//current projection
						Vector2f error=listMatches[i][k].pos_c-x_c;//error in pixels
						Vector2f errorPix=myCamera->m2PixProjJac()*error;//error in pixels
						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
							residue_robust+=TukeyCoef*sqrt(errorPix.squaredNorm());
							valid_wmeas+=TukeyCoef;
								
							miniBAjacz newJc;
							newJc.proj_error=error;
							newJc.de_dp=-TukeyCoef*myCamera->ProjectZ1_Jac_Dp(mapPointsc);
							

							newJc.index_pt_opt=LUT_to_opt[i];
							if(newJc.index_pt_opt!=-1)
							{
								newJc.index_cam_opt=k;
								
								//get jacobien  of error with respect to vairation of pose 2
								Vector3f JacXInvd;
								float inv_depth_cur=idepthp;
								//JacXInvd[0]=-listMatches[i][k].pos_p[0]/(inv_depth_cur*inv_depth_cur);
								//JacXInvd[1]=-listMatches[i][k].pos_p[1]/(inv_depth_cur*inv_depth_cur);
								JacXInvd[0]=-MeanFeatPos[i][0]/(inv_depth_cur*inv_depth_cur);
								JacXInvd[1]=-MeanFeatPos[i][1]/(inv_depth_cur*inv_depth_cur);
								JacXInvd[2]=-1./(inv_depth_cur*inv_depth_cur);
								//JacXInvd=-JacXInvd;
								
								newJc.de_dz=-TukeyCoef*myCamera->ProjectZ1_Jac_X(mapPointsc)*(poseMiniBA[k]).get_rotation()*JacXInvd;							
														
								//if(newJc.de_dz.transpose()*newJc.de_dz>1e-5)//check if point is actually in good configuration
								//{
								//	LUToptToGlob.push_back(id_point);
								//	nbPointsToUpdate++;
								//}
								//else
								//	newJc.index_pt_opt=-1;

								
								fJacobian.push_back(newJc);
								
							}
						}
					}
				}
			}
	
		}
		std::cout<<"residue_robust = "<<residue_robust<<std::endl;
		std::cout<<"valid_wmeas = "<<valid_wmeas<<std::endl;
		
		//compute Hessain and update
		int nbCamsToUpdate=nb_view_saved;
		int nbPointsToUpdate=opt_cpt;
		VectorXf Jtex(nbPointsToUpdate);Jtex.setZero();
		VectorXf Jtep(6*nbCamsToUpdate);Jtep.setZero();

		float Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i]=0;
		float HxxInv[nbPointsToUpdate];
		MatrixXf Hpp(6*nbCamsToUpdate,6*nbCamsToUpdate);Hpp.setZero();
		MatrixXf Hxp(nbPointsToUpdate,6*nbCamsToUpdate);Hxp.setZero();
		
		//get information matrix
		Matrix2f informationMat;
		Matrix2f CovMatrixPix;CovMatrixPix.setZero();CovMatrixPix(0,0)=0.5;CovMatrixPix(1,1)=0.5;
		Matrix2f CovMatrixMet;CovMatrixMet=myCamera->Pix2mProjJac()*CovMatrixPix;
		informationMat.setZero();informationMat(0,0)=1./CovMatrixMet(0,0);informationMat(1,1)=1./CovMatrixMet(1,1);
		
		for(int i=0;i<fJacobian.size();i++)			
		{
			short &pt_opt_id=fJacobian[i].index_pt_opt;
			short &pt_cam_id=fJacobian[i].index_cam_opt;

			//update Jte
			float updatex=fJacobian[i].de_dz.transpose()*informationMat*fJacobian[i].proj_error;
			Jtex[pt_opt_id]-=updatex;
			//update Hessian
			Hxx[pt_opt_id]+=fJacobian[i].de_dz.transpose()*informationMat*fJacobian[i].de_dz;
			
			//update Jte
			VectorXf updatep=fJacobian[i].de_dp.transpose()*informationMat*fJacobian[i].proj_error;
			for(int k=0;k<6;k++)	Jtep[6*pt_cam_id+k]-=updatep[k];
			//update Hessian
			Hpp.block(6*pt_cam_id,6*pt_cam_id,6,6)+=fJacobian[i].de_dp.transpose()*informationMat*fJacobian[i].de_dp;
			
			MatrixXf updatehxp=fJacobian[i].de_dz.transpose()*informationMat*fJacobian[i].de_dp;
			for(int k2=0;k2<6;k2++)Hxp(pt_opt_id,6*pt_cam_id+k2)+=updatehxp(0,k2);						
		}	
		
		for(int i=0;i<nbPointsToUpdate;i++)
			Hxx[i]=(1.+mLMLambda)*Hxx[i];	
		for(int i=0;i<6*nbCamsToUpdate;i++)
			Hpp(i,i)=(1.+mLMLambda)*Hpp(i,i);	
			
		std::vector<int> pt_gone_wrong;
		for(int i=0;i<nbPointsToUpdate;i++)
		{
			HxxInv[i]=1./(Hxx[i]);
		
			if(isnan(HxxInv[i]))
			{
				coutRed<<"Pt gone wring mini BA "<<i<<endlRed;		
				pt_gone_wrong.push_back(LUT_to_main[i]);
			}
		}	
		
		if(pt_gone_wrong.size()!=0)
			coutErrMapping<<"\tBAR: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<pt_gone_wrong.size()<<endlErrMapping;		
		
		MatrixXf HxxInvHxp(nbPointsToUpdate,6*nbCamsToUpdate);//HxxInvHxp=HxxInv*Hxp;
		for(int i=0;i<nbPointsToUpdate;i++)
			for(int c=0;c<nbCamsToUpdate;c++)HxxInvHxp.block(i,6*c,1,6)=HxxInv[i]*Hxp.block(i,6*c,1,6);
		
		VectorXf HxxInvJtp(nbPointsToUpdate);//HxxInvJtp=HxxInv*Jtex;
		for(int i=0;i<nbPointsToUpdate;i++)HxxInvJtp.segment(i,1)=HxxInv[i]*Jtex.segment(i,1);
		
		MatrixXf A(6*nbCamsToUpdate,6*nbCamsToUpdate);	A=Hpp-Hxp.transpose()*HxxInvHxp;
		VectorXf B(6*nbCamsToUpdate);		B=Jtep-Hxp.transpose()*HxxInvJtp;

		Eigen::FullPivLU<MatrixXf> luAc(A);
		VectorXf Dc=luAc.inverse()*B;	

		std::cout<<"apply updates"<<std::endl;
		float gain=0.5;
		//update camera pose
		
		HomogeneousMatrix newposeMiniBA[max_view_in_mini_BA];
		float newinvDepthFeat[nb_feat_ref];
		for(int k=0;k<nb_feat_ref;k++)newinvDepthFeat[k]=invDepthFeat[k];
		
		
		for(int k=0;k<nb_view_saved;k++)
		{
			VectorXf Dci=gain*Dc.segment(6*k,6);
			newposeMiniBA[k]=HomogeneousMatrix(Dci)*poseMiniBA[k];
		}
		
		//update point coord and depths
		VectorXf Dz_all(nbPointsToUpdate);
		Dz_all=HxxInvJtp-HxxInvHxp*Dc;
		bool goneWrong=false;
		for(int id_opt=0;id_opt<nbPointsToUpdate;id_opt++)
		{
			short id_glob=LUT_to_main[id_opt];
			float Dz=gain*Dz_all[id_opt];
			newinvDepthFeat[id_glob]=invDepthFeat[id_glob]+Dz;
			if(newinvDepthFeat[id_glob]<1e-6)newinvDepthFeat[id_glob]=1e-6;
		}
		
		float residue_robustAfter=0;
		float valid_wmeasAfter=0;
		int test=0;
		for(int k=0;k<nb_view_saved;k++)
		{		
			for(int i=0;i<nb_feat_ref;i++)
				if(listMatches[i][k].matched)
			{
				//get corresponding 3d point
				float idepthp=newinvDepthFeat[i];
				if(idepthp!=-1)//depth not initialised yet => use current match to initialise it
				{
					//Vector3f mapPointsc=poseMiniBA[k]*(toHomogeneous(listMatches[i][k].pos_p)/idepthp);
					Vector3f mapPointsc=newposeMiniBA[k]*(toHomogeneous(MeanFeatPos[i])/idepthp);
					if(mapPointsc[2]>0)
					{
						test++;
						Vector2f x_c=myCamera->ProjectZ1(mapPointsc);//current projection
						Vector2f error=listMatches[i][k].pos_c-x_c;//error in pixels
						Vector2f errorPix=myCamera->m2PixProjJac()*error;//error in pixels
						float TukeyCoef=squareRootTukey(errorPix.transpose()*errorPix,sigma_tukey[k]);
						if(TukeyCoef>0)
						{
							residue_robustAfter+=TukeyCoef*sqrt(errorPix.squaredNorm());
							valid_wmeasAfter+=TukeyCoef;
						}
					}
				}
			}
		}
		std::cout<<"test = \t\t\t"<<test<<std::endl;
		std::cout<<"residue_robust/valid_wmeas = \t\t\t"<<residue_robust/valid_wmeas<<std::endl;
		std::cout<<"residue_robustAfter/valid_wmeasAfter = \t"<<residue_robustAfter/valid_wmeasAfter<<std::endl;
		std::cout<<"residue_robustAfter = \t"<<residue_robustAfter<<std::endl;
		std::cout<<"valid_wmeasAfter = \t"<<valid_wmeasAfter<<std::endl;
		
		if(valid_wmeasAfter>0 && residue_robustAfter/valid_wmeasAfter<residue_robust/valid_wmeas)
		{
			for(int k=0;k<nb_feat_ref;k++)invDepthFeat[k]=newinvDepthFeat[k];
			for(int i=0;i<nb_feat_ref;i++)
			{
				//get corresponding 3d point
				float idepthp=invDepthFeat[i];
				int idFeat=indexCandidateFeatureFromVisoId(i);
				if(idFeat!=-1)
					best_features[idFeat].depthInRef=1./idepthp;
			}
			
			
			for(int k=0;k<nb_view_saved;k++)poseMiniBA[k]=newposeMiniBA[k];
			
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
		}
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}
		
	}	
}
*/

inline float det3(float &m00,float &m01,float &m10,float &m11)
{
	return m00*m11-m01*m10;
}
Matrix3f invMatrix3fl(Matrix3f &M)
{
	float invDet=1./M.determinant();
	
	Matrix3f res;
	res(0,0)= det3(M(1,1),M(1,2),M(2,1),M(2,2));res(0,1)=-det3(M(1,0),M(1,2),M(2,0),M(2,2));res(0,2)= det3(M(1,0),M(1,1),M(2,0),M(2,1));
	res(1,0)=-det3(M(0,1),M(0,2),M(2,1),M(2,2));res(1,1)= det3(M(0,0),M(0,2),M(2,0),M(2,2));res(1,2)=-det3(M(0,0),M(0,1),M(2,0),M(2,1));
	res(2,0)= det3(M(0,1),M(0,2),M(1,1),M(1,2));res(2,1)=-det3(M(0,0),M(0,2),M(1,0),M(1,2));res(2,2)= det3(M(0,0),M(0,1),M(1,0),M(1,1));
	return invDet*res;
	
}
void KeyFrame::doMiniBA(std::vector<p_match> &matches,std::vector<uptoscaleFeature> &feats, Camera *myCamera,HomogeneousMatrix22 &pose,bool useLM)
{
	float norm0=sqrt(pose.get_translation().squaredNorm());
	//algebraic solution gives good init
	//matches and feats should be same size and parallel
	int nb_iter=10;
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;
	float gain=0.5;		

	//std::cout<<"matches.size() : "<<matches.size() <<std::endl;
	//std::cout<<"feats.size() : "<<feats.size() <<std::endl;
	
	for(int iter=0;iter<nb_iter;iter++)
	{
		//get Tukey factor
		std::vector<float> vdErrorSquared[2];
		for(int m=0;m<feats.size();m++)
		{
			Vector3f mapPointsCam=feats[m].getLocalCoordinates();
			
			//error in ref
			Vector2f x_d=Vector2f(matches[m].u1p,matches[m].v1p);//desired projection= measurement
			Vector2f x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			Vector2f error=x_d-x_c;
			vdErrorSquared[0].push_back(error.squaredNorm());
			
			//error in cur
			mapPointsCam=pose*mapPointsCam;
			x_d=Vector2f(matches[m].u1c,matches[m].v1c);//desired projection= measurement
			x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			error=x_d-x_c;
			vdErrorSquared[1].push_back(error.squaredNorm());
		}
		float sigma_tukey[2];
		for(int c=0;c<2;c++)
		{
			if(vdErrorSquared[c].size()==0)sigma_tukey[c]=0;
			else
			{
				sigma_tukey[c]=getSigmaSquared(vdErrorSquared[c]);
				if(sigma_tukey[c] < MINSIGMATUKEY)
					sigma_tukey[c] = MINSIGMATUKEY;
			}
		}
		
		//get jacobians
		std::vector<kf_vis_jacobian> fJacobian;
		float residue_robust=0;
		float valid_wmeas=0;
		for(int m=0;m<feats.size();m++)
		{
			Vector3f mapPointsCam=feats[m].getLocalCoordinates();
			
			//error in ref
			Vector2f x_d=Vector2f(matches[m].u1p,matches[m].v1p);//desired projection= measurement
			Vector2f x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			Vector2f error=x_d-x_c;
			float TukeyCoef=squareRootTukey(error.squaredNorm(),sigma_tukey[0]);
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=myCamera->Pix2mProjJac()*(x_d-x_c);
				jac.weight=TukeyCoef*feats[m].recAngle;
				jac.pt_index=m;
				//get jacobian of error with respect to point position
				jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=-1;
				//jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue_robust+=feats[m].recAngle*TukeyCoef*sqrt(error.squaredNorm());
				valid_wmeas+=feats[m].recAngle*TukeyCoef;
			}
			
			//error in cur
			mapPointsCam=pose*mapPointsCam;
			x_d=Vector2f(matches[m].u1c,matches[m].v1c);//desired projection= measurement
			x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			error=x_d-x_c;
			TukeyCoef=squareRootTukey(error.squaredNorm(),sigma_tukey[1]);
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=myCamera->Pix2mProjJac()*(x_d-x_c);
				jac.weight=TukeyCoef*feats[m].recAngle;
				jac.pt_index=m;
				//get jacobian of error with respect to point position
				jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=0;
				jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue_robust+=feats[m].recAngle*TukeyCoef*sqrt(error.squaredNorm());
				valid_wmeas+=feats[m].recAngle*TukeyCoef;
			}
			
		}
		
		//std::cout<<"residue_robust : "<<residue_robust <<std::endl;

		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
		int nbCamsToUpdate=1;
		int nbPointsToUpdate=feats.size();
		VectorXf Jtex(3*nbPointsToUpdate);Jtex.setZero();
		VectorXf Jtep(6*nbCamsToUpdate);Jtep.setZero();

		Matrix3f Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i].setZero();
		Matrix3f HxxInv[nbPointsToUpdate];
		MatrixXf Hpp(6*nbCamsToUpdate,6*nbCamsToUpdate);Hpp.setZero();
		MatrixXf Hxp(3*nbPointsToUpdate,6*nbCamsToUpdate);Hxp.setZero();

		for(int i=0;i<fJacobian.size();i++)			
		{
			int &pt_opt_id=fJacobian[i].pt_index;
			int &cam_opt_id=fJacobian[i].cam_index;

			if(pt_opt_id!=-1)
			{
				//update Jte
				Vector3f updatex=fJacobian[i].weight*fJacobian[i].de_dx.transpose()*fJacobian[i].proj_error;
				for(int k=0;k<3;k++)	Jtex[3*pt_opt_id+k]-=updatex[k];
				//update Hessian
				Hxx[pt_opt_id]+=fJacobian[i].weight*fJacobian[i].de_dx.transpose()*fJacobian[i].de_dx;
			}
			if(cam_opt_id!=-1)
			{
				//update Jte
				VectorXf updatep=fJacobian[i].weight*fJacobian[i].de_dp.transpose()*fJacobian[i].proj_error;
				for(int k=0;k<6;k++)	Jtep[k+6*cam_opt_id]-=updatep[k];
				//update Hessian
				Hpp.block(6*cam_opt_id,6*cam_opt_id,6,6)+=fJacobian[i].weight*fJacobian[i].de_dp.transpose()*fJacobian[i].de_dp;
				
				if(pt_opt_id!=-1)
				{
					MatrixXf updatehxp=fJacobian[i].weight*fJacobian[i].de_dx.transpose()*fJacobian[i].de_dp;
					for(int k=0;k<3;k++)
						for(int k2=0;k2<6;k2++)Hxp(3*pt_opt_id+k,k2+6*cam_opt_id)+=updatehxp(k,k2);						
				}
			}
		}
		/*std::cout<<"Hpp = "<<std::endl;
		std::cout<<Hpp<<std::endl;
		std::cout<<"Hxx = "<<std::endl;
		for(int i=0;i<nbPointsToUpdate;i++)	
			std::cout<<Hxx[i]<<std::endl;*/
		for(int i=0;i<feats.size();i++)
		{
			for(int j=0;j<3;j++)Hxx[i](j,j)=(1.+mLMLambda)*Hxx[i](j,j);
			HxxInv[i]=invMatrix3fl(Hxx[i]);
		}
		
		for(int i=0;i<nbCamsToUpdate;i++)
			for(int j=0;j<6;j++)Hpp(6*i+j,6*i+j)=(1.+mLMLambda)*Hpp(6*i+j,6*i+j);

		MatrixXf HxxInvHxp(3*nbPointsToUpdate,6*nbCamsToUpdate);//HxxInvHxp=HxxInv*Hxp;
		for(int i=0;i<nbPointsToUpdate;i++)
			for(int c=0;c<nbCamsToUpdate;c++)HxxInvHxp.block(3*i,6*c,3,6)=HxxInv[i]*Hxp.block(3*i,6*c,3,6);
		
		VectorXf HxxInvJtp(3*nbPointsToUpdate);//HxxInvJtp=HxxInv*Jtex;
		for(int i=0;i<nbPointsToUpdate;i++)HxxInvJtp.segment(3*i,3)=HxxInv[i]*Jtex.segment(3*i,3);
		
		MatrixXf A(6*nbCamsToUpdate,6*nbCamsToUpdate);	A=Hpp-Hxp.transpose()*HxxInvHxp;
		VectorXf B(6*nbCamsToUpdate);		B=Jtep-Hxp.transpose()*HxxInvJtp;

		Eigen::FullPivLU<MatrixXf> luAc(A);
		VectorXf Dc=luAc.inverse()*B;
		/*std::cout<<"HxxInvHxp = "<<std::endl;
		std::cout<<HxxInvHxp<<std::endl;
		std::cout<<"Hxp = "<<std::endl;
		std::cout<<Hxp<<std::endl;
		std::cout<<"A = "<<std::endl;
		std::cout<<A<<std::endl;*/
		//update camera poses in temporary struct New
		VectorXf Dci=gain*Dc;
        HomogeneousMatrix22 poseNew=HomogeneousMatrix22(Dci)*pose;
		//std::cout<<"Dci = "<<std::endl;
		//std::cout<<Dci.transpose()<<std::endl;
		if(isnan(Dci[0]))
		{
			std::cout<<"nanery in miniBA"<<std::endl;
			break;
		}
		
		//update point positions
		VectorXf Dx_all(3*nbPointsToUpdate);
		Dx_all=HxxInvJtp-HxxInvHxp*Dc;
		//std::cout<<"Dx_all = "<<std::endl;
		//std::cout<<Dx_all<<std::endl;

		Vector3f featsNew[feats.size()];
		for(int i=0;i<feats.size();i++)
		{
			Vector3f Dx;for(int k=0;k<3;k++)Dx[k]=Dx_all[3*i+k];
			featsNew[i]=feats[i].getLocalCoordinates()+gain*Dx;
		}
		
		//compute new residue
		float residue_robust_after=0;
		float valid_wmeas_after=0;
		for(int m=0;m<feats.size();m++)
		{
			Vector3f mapPointsCam=featsNew[m];
			
			//error in ref
			Vector2f x_d=Vector2f(matches[m].u1p,matches[m].v1p);//desired projection= measurement
			Vector2f x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			Vector2f error=x_d-x_c;
			float TukeyCoef=squareRootTukey(error.squaredNorm(),sigma_tukey[0]);
			{
				//update residue to do levendberg marquardt check
				residue_robust_after+=feats[m].recAngle*TukeyCoef*sqrt(error.squaredNorm());
				valid_wmeas_after+=feats[m].recAngle*TukeyCoef;
			}
			
			//error in cur
			mapPointsCam=poseNew*mapPointsCam;
			x_d=Vector2f(matches[m].u1c,matches[m].v1c);//desired projection= measurement
			x_c=myCamera->ToPixels(myCamera->ProjectZ1(mapPointsCam));//current projection				
			error=x_d-x_c;
			TukeyCoef=squareRootTukey(error.squaredNorm(),sigma_tukey[1]);
			{
				//update residue to do levendberg marquardt check
				residue_robust_after+=feats[m].recAngle*TukeyCoef*sqrt(error.squaredNorm());
				valid_wmeas_after+=feats[m].recAngle*TukeyCoef;
			}
			
		}		
		//coutMapping<<"residue before : "<<residue_robust<<"/"<<valid_wmeas<<" = "<<residue_robust/valid_wmeas<<"\t residue after : "<<residue_robust_after<<"/"<<valid_wmeas_after<<" = "<<residue_robust_after/valid_wmeas_after<<endlMapping;

		if(!useLM || residue_robust_after/valid_wmeas_after<residue_robust/valid_wmeas)
		{
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
						
			pose=poseNew;
			for(int i=0;i<feats.size();i++)
			{
				float depth1,depth2;
				depth1=featsNew[i][2];
				depth2=(poseNew*featsNew[i])[2];
				if(depth1<0 || depth2<0)
				{
					matches.erase(matches.begin()+i);
					feats.erase(feats.begin()+i);
					i--;
				}
				else
				{
					feats[i].posRef=myCamera->ProjectZ1(featsNew[i]);
					feats[i].depthInRef=depth1;
				}
			}			
			
		}
		//if not do not confirm change but change lanbda
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}

	}
	
	float normNoScale=sqrt(pose.get_translation().squaredNorm());
	pose.set_translation((norm0/normNoScale)*pose.get_translation());
	for(int i=0;i<feats.size();i++)
	{
		feats[i].depthInRef=(norm0/normNoScale)*feats[i].depthInRef;
	}			
	

}

