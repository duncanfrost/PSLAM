#include "BundleAdjuster.h"

#define MIN_NB_MEASURE_PER_OPTIMISED_CAM 5


inline float det2(float &m00,float &m01,float &m10,float &m11)
{
	return m00*m11-m01*m10;
}

Matrix3f invMatrix3f(Matrix3f &M)
{
	float invDet=1./M.determinant();
	
	Matrix3f res;
	res(0,0)= det2(M(1,1),M(1,2),M(2,1),M(2,2));res(0,1)=-det2(M(1,0),M(1,2),M(2,0),M(2,2));res(0,2)= det2(M(1,0),M(1,1),M(2,0),M(2,1));
	res(1,0)=-det2(M(0,1),M(0,2),M(2,1),M(2,2));res(1,1)= det2(M(0,0),M(0,2),M(2,0),M(2,2));res(1,2)=-det2(M(0,0),M(0,1),M(2,0),M(2,1));
	res(2,0)= det2(M(0,1),M(0,2),M(1,1),M(1,2));res(2,1)=-det2(M(0,0),M(0,2),M(1,0),M(1,2));res(2,2)= det2(M(0,0),M(0,1),M(1,0),M(1,1));
	return invDet*res;
	
}

BundleAdjuster::BundleAdjuster(obMap *_myMap)
{
	myMap=_myMap;
	myCamera=myMap->getCamera();
	nbCamsToUpdate=0;
	nbPointsToUpdate=0;
	hasConverged=false;
	point_translation_small=1e-5;
	cam_translation_small=1e-5;
	goneWrong=false;
#ifdef USE_OMP_C	
	canBeInterrupted=false;
#endif	
	gain=DEFAULT_GAIN_GN;
}
#ifdef USE_OMP_C	
void BundleAdjuster::setMoreImportantStuffToDo(bool *_moreImportantStuffWaiting,omp_lock_t *_lock_check_more_prior)
{
	canBeInterrupted=true;
	moreImportantStuffWaiting=_moreImportantStuffWaiting;
	lock_check_more_prior=_lock_check_more_prior;
}
#endif	

void BundleAdjuster::optimiseInnerWindow(std::vector<int> &innerWin,int nb_iter,bool robust)
{
	std::cout<<"bundleAdjuster optimiseInnerWindow"<<std::endl;
	std::vector<int> optimKF;
	optimKF=innerWin;
	//get outer fix windows:
	std::vector<int> fixedKF=myMap->getDirectNeigbors(optimKF);

	if(optimKF.size()!=0 && fixedKF.size()==0)
	{
		fixedKF.push_back(*optimKF.begin());
		optimKF.erase(optimKF.begin());
	}	
	
	if(optimKF.size()!=0 || fixedKF.size()!=0)
	{
	//create BA:
	
	/*std::cout<<"fixedKF = "<<std::endl;
	for(int i=0;i<fixedKF.size();i++)std::cout<<fixedKF[i]<<"  ";
	std::cout<<std::endl<<std::endl;
	std::cout<<"optimKF = "<<std::endl;
	for(int i=0;i<optimKF.size();i++)std::cout<<optimKF[i]<<"  ";
	std::cout<<std::endl<<std::endl;*/
	
	//pass all the camera positions
	//std::cout<<"add KF"<<std::endl;
	for(int i=0;i<fixedKF.size();i++)
		addCamera(myMap->getKF(fixedKF[i])->getPose(),fixedKF[i],true);
	for(int i=0;i<optimKF.size();i++)
		addCamera(myMap->getKF(optimKF[i])->getPose(),optimKF[i],false);
	
	//need also to pass local best stereo pair of each KF
	for(int i=0;i<fixedKF.size();i++)
		addCamera(myMap->getKF(fixedKF[i])->getBestRelPose()*myMap->getKF(fixedKF[i])->getPose(),-fixedKF[i]-1,false);
	for(int i=0;i<optimKF.size();i++)
		addCamera(myMap->getKF(optimKF[i])->getBestRelPose()*myMap->getKF(optimKF[i])->getPose(),-optimKF[i]-1,false);

	//get points in bundle
	//std::cout<<"add points"<<std::endl;
	std::vector<int> allKF=optimKF;
	for(int i=0;i<fixedKF.size();i++)allKF.push_back(fixedKF[i]);
	
	//have to put all map points (include all matched local features)
	for(int i=0;i<allKF.size();i++)
	{
		KeyFrame &kf=*myMap->getKF(allKF[i]);
		for(int p=0;p<kf.getNbMapPoint();p++)
			if(kf.getPtMapPoint(p)->isUsed())
		{
			addPoint(kf.getPtMapPoint(p)->getPosition(),allKF[i],p,kf.getPtMapPoint(p)->getWeight());
			//std::cout<<"add point: kf :"<<allKF[i]<<" id "<<p<<std::endl;
		}
	}
	//and all local features unmatched
	for(int i=0;i<allKF.size();i++)
	{
		KeyFrame &kf=*myMap->getKF(allKF[i]);
		for(int p=0;p<kf.getNbLocalBestFeatures();p++)
		{
			uptoscaleFeature &feat=*kf.getPtLocalBestFeatures(p);
			if(!feat.matched)
			{
				addPoint(kf.getPose().inverse()*feat.getLocalCoordinates(),allKF[i],-p-1,feat.scoreFundamentalOrigin);
				//std::cout<<"add point: kf :"<<allKF[i]<<" id "<<-p-1<<std::endl;
				
			}
		}
	}
	
	//add measures
	//get measures in Bundle
	//std::cout<<"add measures"<<std::endl;
	for(int c=0;c<allKF.size();c++)
	{
		KeyFrame &kf=*myMap->getKF(allKF[c]);
		//add all measures due to local mini stereo
		std::vector<p_match> &bestMatches=*kf.getBestLocalMatches();
		//std::cout<<"add features from KF "<<allKF[c]<<" bestMatches.size() = "<<bestMatches.size()<<std::endl;
		//at first could think that we could use for loop on LocalBestFeatures but actually will need measure of feature in kf and best stereo=> best matches are better
		for(int m=0;m<bestMatches.size();m++)
		{
			//get corresponding point:
			int id_feat=kf.indexCandidateFeatureFromVisoId(bestMatches[m].i1p);
			if(id_feat==-1)continue;
			uptoscaleFeature &feat=*kf.getPtLocalBestFeatures(id_feat);
		
			int idBA_pt;
			int idBA_kf_pt;
			if(feat.matched)
			{
				//if feature matched then need to link measure to matched point 
				if(!feat.ptKForigin->getPtMapPoint(feat.idPoint)->isUsed())continue;
				idBA_pt=feat.idPoint;
				idBA_kf_pt=feat.ptKForigin-myMap->getKF(0);
			}
			else
			{
				//if not matched then measure is linked to local feature
				idBA_pt=-id_feat-1;
				idBA_kf_pt=allKF[c];
			}
			
			
			//add measure from kf
			//std::cout<<"add meas in "<<allKF[c]<<" and "<<-allKF[c]-1<<std::endl;
			Vector2f measMeterp = myCamera->ToMeters(Vector2f(bestMatches[m].u1p,bestMatches[m].v1p));
			//std::cout<<"\t kf meas : "<<allKF[c]<<" pt: kf "<<idBA_kf_pt<<" id "<<idBA_pt<<std::endl;
			addMeasure(allKF[c],idBA_kf_pt,idBA_pt,measMeterp,0,m);
			
			//add measure from ministereo
			Vector2f measMeterc = myCamera->ToMeters(Vector2f(bestMatches[m].u1c,bestMatches[m].v1c));
			//std::cout<<"\t kf meas : "<<-allKF[c]-1<<" pt: kf "<<idBA_kf_pt<<" id "<<idBA_pt<<std::endl;
			addMeasure(-allKF[c]-1,idBA_kf_pt,idBA_pt,measMeterc,0,m);
			
		}
		
	}
	
	//do bundle adjustment
	//std::cout<<"do BA"<<std::endl;
	optimise(nb_iter,robust);
	
	//retrieve result
	if(getConvergenceResult()!=Diverged)
	{
		std::cout<<"BA converged => update poses"<<std::endl;
		for(int i=0;i<optimKF.size();i++)
		{
			//std::cout<<"before  = "<<myMap->getKF(optimKF[i])->getPose()<<std::endl;
			//std::cout<<"udpated = "<<getUpdatedCamPose(optimKF[i])<<std::endl;
			myMap->getKF(optimKF[i])->setPose(getUpdatedCamPose(optimKF[i]));
		}
		for(int i=0;i<allKF.size();i++)
			myMap->getKF(allKF[i])->setBestRelPose(getUpdatedCamPose(-allKF[i]-1)*myMap->getKF(allKF[i])->getPose().inverse());
		
		
		for(int i=0;i<allKF.size();i++)
		{
			KeyFrame &kf=*myMap->getKF(allKF[i]);
			for(int p=0;p<kf.getNbMapPoint();p++)
				if(kf.getPtMapPoint(p)->isUsed())
					kf.getPtMapPoint(p)->updatePosition(getUpdatedPointPosition(allKF[i],p));
		}
		//and all local features unmatched
		for(int i=0;i<allKF.size();i++)
		{
			KeyFrame &kf=*myMap->getKF(allKF[i]);
			for(int p=0;p<kf.getNbLocalBestFeatures();p++)
			{
				uptoscaleFeature &feat=*kf.getPtLocalBestFeatures(p);
				if(!feat.matched)
				{
					Vector3f newPos=getUpdatedPointPosition(allKF[i],-p-1);
					//project in kf:
					Vector3f localPos=kf.getPose()*newPos;
					feat.depthInRef=localPos[2];
				}
				else
				{
					//get depth from map point
					Vector3f newPos=feat.ptKForigin->getPtMapPoint(feat.idPoint)->getPosition();
					//project in kf:
					Vector3f localPos=kf.getPose()*newPos;
					feat.depthInRef=localPos[2];					
				}
			}
		}
		
		//have to reset relative poses		
		for(int c=0;c<allKF.size();c++)
		{
			KeyFrame &kf=*myMap->getKF(allKF[c]);
			
			for(int n=0;n<kf.getNbNeigbours();n++)
			{
				NeigbourKFNew &neigbor=*kf.getPtNeigbour(n);
				if(std::find(allKF.begin(),allKF.end(),neigbor.neighboring_kf)!=allKF.end())
				{
					//get matches:
					KeyFrame &kfn=*myMap->getKF(neigbor.neighboring_kf);
					HomogeneousMatrix relPose=kf.getPose()*kfn.getPose().inverse();
					neigbor.relative_poses=relPose;
				}
			}
		}
		
		//give outliers to map
		for(int o=0;o<nbOultiers();o++)
		{
			//std::cout<<"remove outlier "<<o<<std::endl;
			OutlierMeasure &mOutlier=getOultiers(o);
			bool outlierOnBestPair=(mOutlier.kf_id_main<0);
			if(outlierOnBestPair)mOutlier.kf_id_main=-(mOutlier.kf_id_main+1);
			
			KeyFrame &kfc=*myMap->getKF(mOutlier.kf_id_main);
			int m=mOutlier.id_feat_main;
			
			p_match &outlier_match=*(kfc.getBestLocalMatches()->begin()+m);
			//get corresponding point:
			int id_feat=kfc.indexCandidateFeatureFromVisoId(outlier_match.i1p);
			if(id_feat==-1)continue;
			uptoscaleFeature &feat=*kfc.getPtLocalBestFeatures(id_feat);
			
			if(feat.matched)
			{
				MapPoint &point=*feat.ptKForigin->getPtMapPoint(feat.idPoint);
				//have to remove view from point 
				point.removeView(kfc.getId());
				
				//if point is now seen by only one other view then need to remove
				//it since does not respect our definition of a point
				if(point.nbViews()==1)
				{
					//unmatched local feature linked to it
					int kfView=point.getView(0);
					int i1pView=point.getI1p(0);
					int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
					if(id_feat_v!=-1)
					{
						uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
						feat_v.matched=false;
					}
					//remove last view
					point.removeViews();
				  
					//remove point
					point.setAsBad();
				}
				
			}
			//remove local feature 
			kfc.removeLocalBestFeatures(id_feat);
			
			
		}
	}
	else
		std::cout<<"BA diverged"<<std::endl;	
	}	
}


void BundleAdjuster::addCamera(HomogeneousMatrix _pose,int _id_main,bool _fixed)
{
	//coutMapping<<"Add Camera "<<_id_main<<endlMapping;
	BAcamera newCam;
	newCam.pose=_pose;
	newCam.poseNew=_pose;
	newCam.nb_measures=0;
	newCam.id_main=_id_main;
	newCam.id_optim=-1;
	newCam.fixed=_fixed;
	mCams.push_back(newCam);
	
	//if(!_fixed)nbCamsToUpdate++;
}

void BundleAdjuster::addPoint(Vector3f _position,int _id_kf_pt,int _id_main,float _recAngle)
{
	//coutMapping<<"Add Point "<<_id_main<<endlMapping;
	BApoint newPoint;
	newPoint.position=_position;
	newPoint.positionNew=_position;
	newPoint.confidence=_recAngle;
	//newPoint.confidence=1.;
	newPoint.nb_measures=0;
	newPoint.id_main=_id_main;
	newPoint.id_kf_main=_id_kf_pt;
	newPoint.id_optim=-1;//will wait to have two measures to put point into optiom
	mPoints.push_back(newPoint);
}

void BundleAdjuster::addMeasure(short _idKFMain,short _idKFPointMain,short _idPointMain,Vector2f _coord,char _lvl,short id_feat)
{
	//coutMapping<<endlMapping<<"Add Measure "<<id_feat<<endlMapping;
	BAmeasure newMeas;
	//get corresponding id of kf in Bundle, could be optimised a bit...
	newMeas.kf_id=-1;
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_idKFMain)
		{
			newMeas.kf_id=c;
			break;
		}
	//if(newMeas.kf_id!=-1)coutMapping<<"mCams["<<newMeas.kf_id<<"].nb_measures "<<mCams[newMeas.kf_id].nb_measures<<endlMapping;

	
	//get corresponding id of point in Bundle, could be optimised a bit...
	newMeas.pt_id=-1;
	for(int p=0;p<mPoints.size();p++)
		if(mPoints[p].id_main==_idPointMain && mPoints[p].id_kf_main==_idKFPointMain)
		{
			newMeas.pt_id=p;
			break;
		}
	//if(newMeas.pt_id!=-1)coutMapping<<"mPoints["<<newMeas.pt_id<<"].nb_measures "<<mPoints[newMeas.pt_id].nb_measures<<endlMapping;
	

	//will wait to have two measures to put point into optiom
	if(newMeas.kf_id!=-1 && newMeas.pt_id!=-1)//kf_id ==-1 should cetrainly not happen but well, never too carefull
	{
		mCams[newMeas.kf_id].nb_measures++;
		mPoints[newMeas.pt_id].nb_measures++;
		
		//check of cam is not fixed and had enough measures to be optimised
		if(mCams[newMeas.kf_id].fixed==false && mCams[newMeas.kf_id].nb_measures==MIN_NB_MEASURE_PER_OPTIMISED_CAM)
		{
			mCams[newMeas.kf_id].id_optim=nbCamsToUpdate;
			nbCamsToUpdate++;
		}
	
		//check of cam is not fixed and had enough measures to be optimised
		if(mPoints[newMeas.pt_id].nb_measures==2)
		{
			mPoints[newMeas.pt_id].id_optim=nbPointsToUpdate;
			nbPointsToUpdate++;
		}
		
		newMeas.coord=_coord;
		newMeas.lvl_meas=_lvl;
		newMeas.id_feat_main =id_feat;
		mMeasures.push_back(newMeas);
	}
}


void BundleAdjuster::optimise(int nb_iter,bool robust)
{
	//coutMapping<<"BundleAdjuster::optimise> nbCamsToUpdate : "<<nbCamsToUpdate<<endlMapping;
	//set cam to be fixed if they don t have enough measures
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].nb_measures<MIN_NB_MEASURE_PER_OPTIMISED_CAM)mCams[c].fixed=true;//is not really important there cause we can check id_optim==-1 but still its cleaner like this
	
	//for(int i=0;i<mCams.size();i++)
	//	coutMapping<<"mCams "<<i<<" : nbMeas = "<<mCams[i].nb_measures<<endlMapping;
	//for(int i=0;i<mPoints.size();i++)
	//	coutMapping<<"mPoints "<<i<<" : nbMeas = "<<mPoints[i].nb_measures<<" id_optim = "<<mPoints[i].id_optim<<endlMapping;

	//Check how many cameras have to be updated, if zero => optimise point position only
	if(nbCamsToUpdate==0)
	{
#ifdef VERBOSE
		coutMapping<<"\tBundleAdjuster::optimise> nbCamsToUpdate : "<<nbCamsToUpdate<<" => only do point optim"<<endlMapping;
#endif
        OptimisePointPosition(nb_iter);
	}
	//if more than 1 then can do proper BA (example : one cam fixed and one cam to be updated)
	else
	{
		if(robust)
		{
#ifdef VERBOSE
			coutMapping<<"\tBundleAdjuster::optimise> nbCamsToUpdate : "<<nbCamsToUpdate<<" => do BundleAdjustRobust"<<endlMapping;
#endif
			BundleAdjustRobust2(nb_iter);
		}
		else
		{
#ifdef VERBOSE
			coutMapping<<"\tBundleAdjuster::optimise> nbCamsToUpdate : "<<nbCamsToUpdate<<" => do BundleAdjust"<<endlMapping;
#endif
            BundleAdjust(nb_iter,true);
		}
	}
}
void BundleAdjuster::OptimisePointPosition(int nb_iter)
{
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();

	for(int iter=0;iter<nb_iter;iter++)
	{
		
		//if has done one iteration then position of point should be good enough to be checked
		if(iter>0)
			for(int m=0;m<mMeasures.size();m++)
			{
				short &c=mMeasures[m].kf_id;
				short &p=mMeasures[m].pt_id;
				if(mPoints[p].id_optim!=-1)
				{
					HomogeneousMatrix &pose=mCams[c].pose;
					
					//get 3D pose in camera frame
					Vector3f mapPointsCam=pose*mPoints[p].position;	
					
					//compute error with observed points (in meter in z=1 plane)
					Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
					Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				

					//compute jacobian only if point or cam linked to measure is to be optimised
					Vector2f proj_error=x_d-x_c;
					
					float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*proj_error).squaredNorm());
					float TukeyCoef=squareRootTukey(norm_reproj_error,4);//if its more than 2 px then has to be outlier
					
					if(iter>0 && TukeyCoef==0)//measure is considered as outlier
					{
						rejectMeasure(m);
						m--;
					}
				}
			}
		
		//compute all non null elements of jacobien and reprojection error
		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
		Vector3f Jtex[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Jtex[i].setZero();

		Matrix3f Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i].setZero();
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			if(mPoints[p].id_optim!=-1)
			{
				HomogeneousMatrix &pose=mCams[c].pose;
				
				//get 3D pose in camera frame
				Vector3f mapPointsCam=pose*mPoints[p].position;	
				
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
				Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				

				//compute jacobian only if point or cam linked to measure is to be optimised
				Vector2f proj_error=x_d-x_c;
				
				//get jacobian of error with respect to point position
				MatrixXf de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				int &pt_opt_id=mPoints[p].id_optim;
				//update Jte
				Vector3f updatex=de_dx.transpose()*proj_error;
				for(int k=0;k<3;k++)	Jtex[pt_opt_id][k]-=updatex[k];
				//update Hessian
				Hxx[pt_opt_id]+=de_dx.transpose()*de_dx;			
			}
		}

		int cpt_pt_no_dev=0;
		float max_point_translation_update=0;
		for(int p=0;p<mPoints.size();p++)
		{
			int &pt_opt_id=mPoints[p].id_optim;
			if(pt_opt_id!=-1)
			{
				//Eigen::FullPivLU<MatrixXf> luHxxi(Hxx[i]);
				//Matrix3f HxxInv=luHxxi.inverse();
				//for(int j=0;j<3;j++)Hxx[i](j,j)+=LevMarLikeConstant;
				Matrix3f HxxInv=invMatrix3f(Hxx[pt_opt_id]);

				if(isnan(HxxInv(0,0)))
				{
					cpt_pt_no_dev++;
				}
				else
				{
					Vector3f Dx=gain*HxxInv*Jtex[pt_opt_id];
					
					mPoints[p].positionNew=mPoints[p].position+Dx;
					if(!isnan(mPoints[p].positionNew[0]))//could happen if wrong configuration (point toward infinity or wrong measures)
					{
						mPoints[p].position=mPoints[p].positionNew;
						float translation_update=Dx.squaredNorm();
						if(max_point_translation_update<translation_update)
							max_point_translation_update=translation_update;
					}
					else
					{
						mPoints[p].positionNew=mPoints[p].position;
						//do not optimise point anymore
						//doNotOptimisePoint(p);
						rejectPoint(p);
					}

					
				}	
			}
			
		}
		
		if(cpt_pt_no_dev!=0)
			coutErrMapping<<"\tOptimPointPos: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<cpt_pt_no_dev<<endlErrMapping;
		
		//convergence test
		if(max_point_translation_update<point_translation_small)
		{
			hasConverged=true;
			break;
		}
#ifdef USE_OMP_C		
		if(canBeInterrupted)
		{
			omp_set_lock(lock_check_more_prior);
			if(*moreImportantStuffWaiting)
			{
			  omp_unset_lock(lock_check_more_prior);	
			  break; 
			}			
			omp_unset_lock(lock_check_more_prior);	
		}
#endif	

		
	}//end of for iter
}

void BundleAdjuster::rejectPoint(int _p)
{
	for(int m=0;m<mMeasures.size();m++)
	{
		short &p=mMeasures[m].pt_id;
		if(p==_p)
		{
			rejectMeasure(m);
			m--;
		}
	}
}

void BundleAdjuster::doNotOptimisePoint(int p)
{
	if(mPoints[p].id_optim!=-1)
	{
		//set it as fixed
		//=>all the other points with id_optim > mPoints[p].id_optim have to be shifted
		for(int p2=0;p2<mPoints.size();p2++)
			if(mPoints[p2].id_optim>mPoints[p].id_optim)
				mPoints[p2].id_optim--;
			
		mPoints[p].id_optim=-1;
		nbPointsToUpdate--;		
	}
}


void BundleAdjuster::rejectMeasure(int m)
{
	short &c=mMeasures[m].kf_id;
	short &p=mMeasures[m].pt_id;
	//decrease number of measure cam and point
	mCams[c].nb_measures--;
	mPoints[p].nb_measures--;
	
	//check if cam and point are still differentiable
	if(mCams[c].id_optim!=-1 && mCams[c].nb_measures<MIN_NB_MEASURE_PER_OPTIMISED_CAM)//not differentiable anymore
	{
		//set it as fixed
		//=>all the other cam with id_optim > mCams[c].id_optim have to be shifted
		for(int c2=0;c2<mCams.size();c2++)
			if(mCams[c2].id_optim>mCams[c].id_optim)
				mCams[c2].id_optim--;
			
		mCams[c].id_optim=-1;
		mCams[c].fixed=true;
		nbCamsToUpdate--;
	}
	
	if(mPoints[p].id_optim!=-1 && mPoints[p].nb_measures<2)//not differentiable anymore set it as fixed
		doNotOptimisePoint(p);
	
	
	//add measure to outlier list
	OutlierMeasure newOutlier;
	newOutlier.kf_id_main=mCams[c].id_main;
	newOutlier.lvl_meas=mMeasures[m].lvl_meas;
	newOutlier.id_feat_main=mMeasures[m].id_feat_main;
	mOutliers.push_back(newOutlier);
	
	//remove measure
	mMeasures.erase(mMeasures.begin()+m);	
}

int BundleAdjuster::IdOptimToIdBA(int i)
{
	for(int i2=0;i2<mPoints.size();i2++)
	{
		if(i==mPoints[i2].id_optim)
			return i2;
	}
	return -1; 
}

void BundleAdjuster::BundleAdjustRobust(int nb_iter)
{
	
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();
	
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;

	for(int iter=0;iter<nb_iter;iter++)
	{
		//std::cout<<"iteration : "<<iter <<std::endl;
		startIterRobust:
		
		//get Tukey factor
		std::vector<float> vdErrorSquared;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			if(mapPointsCam[2]>0)
			{
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
				Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
				Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
				vdErrorSquared.push_back(error.transpose()*error);
			}
		}
		float sigma_tukey=getSigmaSquared(vdErrorSquared);
		if(sigma_tukey < MINSIGMATUKEY)
			sigma_tukey = MINSIGMATUKEY;

		
		
		//points will be considered for optimisation only if they have two valid projection
		//that was kind of done already but Tukey function was not known yet
		//=> update taking robust function into account
		//std::cout<<"nb measure before outlier check : "<<mMeasures.size() <<std::endl;
		std::cout<<"nbCamsToUpdate : "<<nbCamsToUpdate <<std::endl;
		//std::cout<<"sigma_tukey : "<<sigma_tukey <<std::endl;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			  
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);
			//std::cout<<"kf :"<<c <<"  error = "<<(myCamera->m2PixProjJac()*(x_d-x_c)).transpose()<<std::endl;
			//std::cout<<"norm_reproj_error :"<<norm_reproj_error <<"  TukeyCoef = "<<TukeyCoef<<std::endl;
			if(mapPointsCam[2]<0 || mPoints[p].confidence*TukeyCoef==0)//measure is considered as outlier
			{
				//std::cout<<"reject measure from kf :"<<c <<std::endl;
				rejectMeasure(m);
				m--;
			}
		}
		
		//std::cout<<"nb measure after outlier check : "<<mMeasures.size() <<std::endl;
		std::cout<<"nbCamsToUpdate after: "<<nbCamsToUpdate <<std::endl;
		
		//check if BA is still possible
		if(nbCamsToUpdate==0)
		{
			OptimisePointPosition(nb_iter-iter);
			break;//has done point optim => goto end
		}
		
		//compute all non null elements of jacobien and reprojection error
		std::vector<kf_vis_jacobian> fJacobian;
		float residue_robust=0;
		float valid_wmeas=0;
		
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);
			
			//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mapPointsCam[2]>0 && mPoints[p].confidence*TukeyCoef!=0 && (mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1))
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=x_d-x_c;
				jac.weight=TukeyCoef*mPoints[p].confidence;
				jac.pt_index=mPoints[p].id_optim;
				//get jacobian of error with respect to point position
				if(mPoints[p].id_optim!=-1)
					jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=mCams[c].id_optim;
				if(mCams[c].id_optim!=-1)
					jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue_robust+=mPoints[p].confidence*TukeyCoef*sqrt(norm_reproj_error);
				valid_wmeas+=mPoints[p].confidence*TukeyCoef;
				
				//if(isnan(jac.proj_error[0]))std::cout<<"nanery in error "<<std::endl;
			}
		}
		std::cout<<"residue_robust : "<<residue_robust <<std::endl;
		//std::cout<<"residue_robust norm: "<<residue_robust/valid_wmeas <<std::endl;
		

		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
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
		
		
		std::vector<int> pt_gone_wrong;
		for(int i=0;i<nbPointsToUpdate;i++)
		{
			for(int j=0;j<3;j++)Hxx[i](j,j)=(1.+mLMLambda)*Hxx[i](j,j);
			//Eigen::FullPivLU<MatrixXf> luHxxi(Hxx[i]);
			//HxxInv[i]=luHxxi.inverse();
			HxxInv[i]=invMatrix3f(Hxx[i]);

			
			if(isnan(HxxInv[i](0,0)))
			{
				//find corresponding point and put in list not to be optimised anymore (if it is point at infinity then could still be used to locate camera pose)
				int id_pt=IdOptimToIdBA(i);
				pt_gone_wrong.push_back(id_pt);
				
			}
		}
		if(pt_gone_wrong.size()!=0)
		{
			coutErrMapping<<"\tBAR: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<pt_gone_wrong.size()<<endlErrMapping;
			for(int i=0;i<pt_gone_wrong.size();i++)
				doNotOptimisePoint(pt_gone_wrong[i]);
			
			//go back to start iteration loop
			goto startIterRobust;
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

		//update camera poses in temporary struct New

		float max_cam_translation_update=0;
		for(int c=0;c<nb_keyFrames;c++)
		{
			short id_opt=mCams[c].id_optim;
			if(id_opt!=-1)
			{
				
				VectorXf Dci=gain*Dc.segment(6*id_opt,6);
				mCams[c].poseNew=HomogeneousMatrix(Dci)*mCams[c].pose;
				
				if(isnan(Dci[0]))
				{
					//coutErrMapping<<"\tBA: strange, cam update"<< iter <<endlErrMapping;
					mCams[c].poseNew=mCams[c].pose;
					if(!goneWrong)
						std::cout<<"Dci["<<id_opt<<"] is nan" <<std::endl;
					goneWrong=true;
				}
				
				float translation_update=(Dci.segment(0,3)).squaredNorm();
				if(max_cam_translation_update<translation_update)
					max_cam_translation_update=translation_update;
			}
		}
		//std::cout<<"max_translation_update = "<<max_translation_update<<std::endl;


		//update point positions
		VectorXf Dx_all(3*nbPointsToUpdate);
		Dx_all=HxxInvJtp-HxxInvHxp*Dc;

		float max_point_translation_update=0;
		for(int i=0;i<nb_points;i++)
		{
			short id_opt=mPoints[i].id_optim;
			if(id_opt!=-1)
			{
				Vector3f Dx;for(int k=0;k<3;k++)Dx[k]=gain*Dx_all[3*id_opt+k];
				mPoints[i].positionNew=mPoints[i].position+Dx;
				
				if(isnan(Dx[0]))
				{
					//coutErrMapping<<"\tBA: strange, point update"<< iter <<endlErrMapping;
					mPoints[i].positionNew=mPoints[i].position;
					if(!goneWrong)
						std::cout<<"Dx["<<id_opt<<"] is nan" <<std::endl;
					goneWrong=true;
				}
				
				float translation_update=Dx.squaredNorm();
				if(max_point_translation_update<translation_update)
					max_point_translation_update=translation_update;
			}
		}
		
		//convergence test
		if(max_point_translation_update<point_translation_small && max_cam_translation_update<cam_translation_small)
		{
			hasConverged=true;
			break;
		}		
				
		//check if step was good
		float residue_robust_after=0;
		float valid_meas_after=0;
		int nb_outliers_new=0;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].poseNew;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].positionNew;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey);
			
			//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mPoints[p].confidence*TukeyCoef!=0 && (mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1))
			{
				//residue to do levendberg marquardt check
				residue_robust_after+=mPoints[p].confidence*TukeyCoef*sqrt(norm_reproj_error);
				valid_meas_after+=mPoints[p].confidence*TukeyCoef;			}
			else
				nb_outliers_new++;
		}	
		
		coutMapping<<"residue before : "<<residue_robust<<"/"<<valid_wmeas<<" = "<<residue_robust/valid_wmeas<<"\t residue after : "<<residue_robust_after<<"/"<<valid_meas_after<<" = "<<residue_robust_after/valid_meas_after<<endlMapping;
					
		
		//if step was good confirm change and change lanbda
		if(residue_robust_after/valid_meas_after<residue_robust/valid_wmeas)
		{
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
			
			
			for(int c=0;c<nb_keyFrames;c++)
				mCams[c].pose=mCams[c].poseNew;
			for(int i=0;i<nb_points;i++)
				mPoints[i].position=mPoints[i].positionNew;
			
			if(nb_outliers_new!=0 && iter==nb_iter-1)
				//std::cout<<"There was "<<nb_outliers_new<<"new outliers"<<std::endl;
				coutColMapping<<"There was "<<nb_outliers_new<<"new outliers"<<endlColMapping;
			
			
		}
		//if not do not confirm change but change lanbda
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}
#ifdef USE_OMP_C			
		if(canBeInterrupted)
		{
			omp_set_lock(lock_check_more_prior);
			if(*moreImportantStuffWaiting)
			{
			  omp_unset_lock(lock_check_more_prior);	
			  break; 
			}
			omp_unset_lock(lock_check_more_prior);	
		}
#endif	

		
	}//end of for iter



	//update all relative poses (dp not need to do it for each iteration)
	//std::vector<int> innerWin;for(int i=0;i<KeyFrameList.size();i++)innerWin.push_back(i);
	//updateRelativePoseInWin(innerWin);
	
}


HomogeneousMatrix BundleAdjuster::getUpdatedCamPose(int id_main)
{
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==id_main)return mCams[c].pose;
		
	std::cerr<<"BundleAdjuster::getUpdatedCamPose> cam not found"<<endlMapping;
	return HomogeneousMatrix();//should not happen
}
Vector3f BundleAdjuster::getUpdatedPointPosition(int id_kf_main,int id_main)
{
	for(int p=0;p<mPoints.size();p++)
		if(mPoints[p].id_main==id_main && mPoints[p].id_kf_main==id_kf_main)return mPoints[p].position;
		
	std::cerr<<"BundleAdjuster::getUpdatedPointPosition> point not found"<<endlMapping;
	Vector3f defaultPoint;
	return defaultPoint;
	
}

void BundleAdjuster::BundleAdjust(int nb_iter,bool LM)
{
	
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();
	
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;

	for(int iter=0;iter<nb_iter;iter++)
	{
		startIter:

		std::cout<<"iteration "<<iter<<std::endl;
		//compute all non null elements of jacobien and reprojection error
		std::vector<kf_vis_jacobian> fJacobian;
		float residue=0;
		//float residue_pix=0;
		float valid_wmeas=0;

		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1)
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=x_d-x_c;
				jac.weight=mPoints[p].confidence;
				jac.pt_index=mPoints[p].id_optim;
				//get jacobian of error with respect to point position
				if(mPoints[p].id_optim!=-1)
					jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=mCams[c].id_optim;
				if(mCams[c].id_optim!=-1)
					jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue+=mPoints[p].confidence*sqrt(norm_reproj_error);
				valid_wmeas++;
				
				//just for debug
				/*float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
				std::cout<<"residue cam ["<<mCams[c].id_main<<"] : "<< norm_reproj_error<<std::endl;
				residue_pix+=norm_reproj_error;*/
				
			}
		}
		//std::cout<<"residue_pix_norm ="<<residue_pix/valid_wmeas<<std::endl;
		//std::cout<<"residue ="<<residue<<std::endl;
		

		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
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
			//coutMapping<<"pt_opt_id = "<<pt_opt_id<<endlMapping;

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
		std::vector<int> pt_gone_wrong;
		for(int i=0;i<nbPointsToUpdate;i++)
		{
			if(LM)
				for(int j=0;j<3;j++)Hxx[i](j,j)=(1.+mLMLambda)*Hxx[i](j,j);
			else
				for(int j=0;j<3;j++)Hxx[i](j,j)+=LevMarLikeConstant;
	
			//Eigen::FullPivLU<MatrixXf> luHxxi(Hxx[i]);
			//HxxInv[i]=luHxxi.inverse();
			HxxInv[i]=invMatrix3f(Hxx[i]);

			if(isnan(HxxInv[i](0,0)))
			{
				int id_pt=IdOptimToIdBA(i);
				pt_gone_wrong.push_back(id_pt);
			}
		}
		if(pt_gone_wrong.size()!=0)
		{
			coutErrMapping<<"\tBA: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<pt_gone_wrong.size()<<endlErrMapping;
			for(int i=0;i<pt_gone_wrong.size();i++)
				doNotOptimisePoint(pt_gone_wrong[i]);
			//go back to start iteration loop
			goto startIter;
		}			
		if(LM)
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

		//update camera poses in temporary struct New

		float max_cam_translation_update=0;
		for(int c=0;c<nb_keyFrames;c++)
		{
			short id_opt=mCams[c].id_optim;
			if(id_opt!=-1)
			{
				
				VectorXf Dci=gain*Dc.segment(6*id_opt,6);
				mCams[c].poseNew=HomogeneousMatrix(Dci)*mCams[c].pose;
				if(isnan(Dci[0]))
				{
					//coutErrMapping<<"\tBA: strange, cam update"<< iter <<endlErrMapping;
					mCams[c].poseNew=mCams[c].pose;
					goneWrong=true;
				}

				
				if(!LM)mCams[c].pose=mCams[c].poseNew;
				
				float translation_update=(Dci.segment(0,3)).squaredNorm();
				if(max_cam_translation_update<translation_update)
					max_cam_translation_update=translation_update;
			}
		}
		//std::cout<<"max_cam_translation_update = "<<max_cam_translation_update<<std::endl;

		//update point positions
		VectorXf Dx_all(3*nbPointsToUpdate);
		Dx_all=HxxInvJtp-HxxInvHxp*Dc;

		float max_point_translation_update=0;
		for(int i=0;i<nb_points;i++)
		{
			short id_opt=mPoints[i].id_optim;
			if(id_opt!=-1)
			{
				Vector3f Dx;for(int k=0;k<3;k++)Dx[k]=gain*Dx_all[3*id_opt+k];
				mPoints[i].positionNew=mPoints[i].position+Dx;	
				if(isnan(Dx[0]))
				{
					//coutErrMapping<<"\tBA: strange, cam update"<< iter <<endlErrMapping;
					mPoints[i].positionNew=mPoints[i].position;
					goneWrong=true;
				}

				
				if(!LM)mPoints[i].position=mPoints[i].positionNew;
				
				float translation_update=Dx.squaredNorm();
				if(max_point_translation_update<translation_update)
					max_point_translation_update=translation_update;

			}
		}
		//std::cout<<"max_point_translation_update = "<<max_point_translation_update<<std::endl;
		//convergence test
		if(max_point_translation_update<point_translation_small && max_cam_translation_update<cam_translation_small)
		{
			hasConverged=true;
			break;
		}
				
		if(LM)
		{
			//check if step was good
			float residue_after=0;
			float valid_meas_after=0;
			for(int m=0;m<mMeasures.size();m++)
			{
				short &c=mMeasures[m].kf_id;
				short &p=mMeasures[m].pt_id;
				HomogeneousMatrix &pose=mCams[c].poseNew;
				
				//get 3D pose in camera frame
				Vector3f mapPointsCam=pose*mPoints[p].positionNew;	
				
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
				Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
				//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
				float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
				
				//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
				//compute jacobian only if point or cam linked to measure is to be optimised
				if(mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1)
				{
					//residue to do levendberg marquardt check
					residue_after+=mPoints[p].confidence*sqrt(norm_reproj_error);
					valid_meas_after++;
				}
			}	
			
			//coutMapping<<"residue before : "<<residue<<"/"<<valid_wmeas<<" = "<<residue/valid_wmeas<<"\t residue after : "<<residue_after<<"/"<<valid_meas_after<<" = "<<residue_after/valid_meas_after<<endlMapping;
						
			
			//if step was good confirm change and change lanbda
			if(residue_after/valid_meas_after<residue/valid_wmeas)
			{
				mdLambdaFactor = 2.0;
				mLMLambda *= 0.3;
				
				
				for(int c=0;c<nb_keyFrames;c++)
					mCams[c].pose=mCams[c].poseNew;
				for(int i=0;i<nb_points;i++)
					mPoints[i].position=mPoints[i].positionNew;
				
				
			}
			//if not do not confirm change but change lanbda
			else
			{
				mLMLambda = mLMLambda * mdLambdaFactor;
				mdLambdaFactor = mdLambdaFactor * 2;
			}
					
		}
#ifdef USE_OMP_C			
		if(canBeInterrupted)
		{
			omp_set_lock(lock_check_more_prior);
			if(*moreImportantStuffWaiting)
			{
			  omp_unset_lock(lock_check_more_prior);	
			  break; 
			}
			omp_unset_lock(lock_check_more_prior);	
		}
#endif	
	}//end of for iter

	
}
void BundleAdjuster::BundleAdjustNoOptim(int nb_iter)
{
	
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();
	
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;

	for(int iter=0;iter<nb_iter;iter++)
	{

		
		//compute all non null elements of jacobien and reprojection error
		std::vector<kf_vis_jacobian> fJacobian;
		float residue=0;
		float valid_wmeas=0;

		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1)
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=x_d-x_c;
				jac.pt_index=mPoints[p].id_optim;
				//get jacobian of error with respect to point position
				if(mPoints[p].id_optim!=-1)
					jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=mCams[c].id_optim;
				if(mCams[c].id_optim!=-1)
					jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue+=sqrt(norm_reproj_error);
				valid_wmeas++;
			}
		}
		

		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
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
			//coutMapping<<"pt_opt_id = "<<pt_opt_id<<endlMapping;

			if(pt_opt_id!=-1)
			{
				//update Jte
				Vector3f updatex=fJacobian[i].de_dx.transpose()*fJacobian[i].proj_error;
				for(int k=0;k<3;k++)	Jtex[3*pt_opt_id+k]-=updatex[k];
				//update Hessian
				Hxx[pt_opt_id]+=fJacobian[i].de_dx.transpose()*fJacobian[i].de_dx;
			}
			if(cam_opt_id!=-1)
			{
				//update Jte
				VectorXf updatep=fJacobian[i].de_dp.transpose()*fJacobian[i].proj_error;
				for(int k=0;k<6;k++)	Jtep[k+6*cam_opt_id]-=updatep[k];
				//update Hessian
				Hpp.block(6*cam_opt_id,6*cam_opt_id,6,6)+=fJacobian[i].de_dp.transpose()*fJacobian[i].de_dp;
				
				if(pt_opt_id!=-1)
				{
					MatrixXf updatehxp=fJacobian[i].de_dx.transpose()*fJacobian[i].de_dp;
					for(int k=0;k<3;k++)
						for(int k2=0;k2<6;k2++)Hxp(3*pt_opt_id+k,k2+6*cam_opt_id)+=updatehxp(k,k2);						
				}
			}
		}
		
		MatrixXf H(3*nbPointsToUpdate+6*nbCamsToUpdate,3*nbPointsToUpdate+6*nbCamsToUpdate);H.setZero();
		VectorXf J(3*nbPointsToUpdate+6*nbCamsToUpdate);
		
		for(int i=0;i<3*nbPointsToUpdate;i++)J[i]=Jtex[i];
		for(int i=0;i<6*nbCamsToUpdate;i++)J[3*nbPointsToUpdate+i]=Jtep[i];
	
		for(int k=0;k<nbPointsToUpdate;k++)
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)H(3*k+i,3*k+j)=Hxx[k](i,j);
				
		for(int i=0;i<6*nbCamsToUpdate;i++)
			for(int j=0;j<6*nbCamsToUpdate;j++)H(3*nbPointsToUpdate+i,3*nbPointsToUpdate+j)=Hpp(i,j);
		for(int i=0;i<3*nbPointsToUpdate;i++)
			for(int j=0;j<6*nbCamsToUpdate;j++)
			{
				H(i,3*nbPointsToUpdate+j)=Hxp(i,j);
				H(3*nbPointsToUpdate+j,i)=Hxp(i,j);
			}
		
		//compute update
		Eigen::FullPivLU<MatrixXf> lu(H);
		VectorXf Dv(3*nbPointsToUpdate+6*nbCamsToUpdate);
		Dv=1.0*(lu.inverse()*J);

		//update cam positions
		VectorXf Dc(6*nbCamsToUpdate);
		Dc=Dv.segment(3*nbPointsToUpdate,3*nbPointsToUpdate+6*nbCamsToUpdate);


		//update camera poses in temporary struct New

		for(int c=0;c<nb_keyFrames;c++)
		{
			short id_opt=mCams[c].id_optim;
			if(id_opt!=-1)
			{
				
				VectorXf Dci=gain*Dc.segment(6*id_opt,6);
				mCams[c].poseNew=HomogeneousMatrix(Dci)*mCams[c].pose;
			}
		}


		//update point positions
		VectorXf Dx_all(3*nbPointsToUpdate);
		Dx_all=Dv.segment(0,3*nbPointsToUpdate);

		for(int i=0;i<nb_points;i++)
		{
			short id_opt=mPoints[i].id_optim;
			if(id_opt!=-1)
			{
				Vector3f Dx;for(int k=0;k<3;k++)Dx[k]=gain*Dx_all[3*id_opt+k];
				mPoints[i].positionNew=mPoints[i].position+Dx;
			}
		}
				
		//check if step was good
		float residue_after=0;
		float valid_meas_after=0;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].poseNew;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].positionNew;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			
			//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1)
			{
				//residue to do levendberg marquardt check
				residue_after+=sqrt(norm_reproj_error);
				valid_meas_after++;
			}
		}	
		
		//coutMapping<<"residue before : "<<residue<<"/"<<valid_wmeas<<" = "<<residue/valid_wmeas<<"\t residue after : "<<residue_after<<"/"<<valid_meas_after<<" = "<<residue_after/valid_meas_after<<endlMapping;
					
		
		//if step was good confirm change and change lanbda
		if(residue_after/valid_meas_after<residue/valid_wmeas)
		{
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
			
			
			for(int c=0;c<nb_keyFrames;c++)
				mCams[c].pose=mCams[c].poseNew;
			for(int i=0;i<nb_points;i++)
				mPoints[i].position=mPoints[i].positionNew;
			
			
		}
		//if not do not confirm change but change lanbda
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}

		
	}//end of for iter
	
}

ConvergenceResult BundleAdjuster::getConvergenceResult()
{
	if(goneWrong)return Diverged;
	if(hasConverged)return Converged;
	else return NotConverged;	
};

/*void BundleAdjuster::OptimisePointPosition(int nb_iter)
{
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();

	for(int iter=0;iter<nb_iter;iter++)
	{
		//compute all non null elements of jacobien and reprojection error
		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
		Vector3f Jtex[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Jtex[i].setZero();

		Matrix3f Hxx[nbPointsToUpdate];for(int i=0;i<nbPointsToUpdate;i++)Hxx[i].setZero();
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			if(mPoints[p].id_optim!=-1)
			{
				HomogeneousMatrix &pose=mCams[c].pose;
				
				//get 3D pose in camera frame
				Vector3f mapPointsCam=pose*mPoints[p].position;	
				
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
				Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				

				//compute jacobian only if point or cam linked to measure is to be optimised
				Vector2f proj_error=x_d-x_c;
				//get jacobian of error with respect to point position
				MatrixXf de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				int &pt_opt_id=mPoints[p].id_optim;
				//update Jte
				Vector3f updatex=de_dx.transpose()*proj_error;
				for(int k=0;k<3;k++)	Jtex[pt_opt_id][k]-=updatex[k];
				//update Hessian
				Hxx[pt_opt_id]+=de_dx.transpose()*de_dx;
			}
		}

		int cpt_pt_no_dev=0;
		float max_point_translation_update=0;
		for(int p=0;p<mPoints.size();p++)
		{
			int &pt_opt_id=mPoints[p].id_optim;
			if(pt_opt_id!=-1)
			{
				//Eigen::FullPivLU<MatrixXf> luHxxi(Hxx[i]);
				//Matrix3f HxxInv=luHxxi.inverse();
				//for(int j=0;j<3;j++)Hxx[i](j,j)+=LevMarLikeConstant;
				Matrix3f HxxInv=invMatrix3f(Hxx[pt_opt_id]);

				if(isnan(HxxInv(0,0)))
				{
					cpt_pt_no_dev++;
				}
				else
				{
					Vector3f Dx=HxxInv*Jtex[pt_opt_id];
					
					mPoints[p].positionNew=mPoints[p].position+Dx;
					if(!isnan(mPoints[p].positionNew[0]))//could happen if wrong configuration (point toward infinity or wrong measures)
					{
						mPoints[p].position=mPoints[p].positionNew;
						float translation_update=Dx.squaredNorm();
						if(max_point_translation_update<translation_update)
							max_point_translation_update=translation_update;
					}
					else
					{
						mPoints[p].positionNew=mPoints[p].position;
						//do not optimise point anymore
						//doNotOptimisePoint(p);
						rejectPoint(p);
					}

					
				}	
			}
			
		}
		
		if(cpt_pt_no_dev!=0)
			coutErrMapping<<"\tOptimPointPos: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<cpt_pt_no_dev<<endlErrMapping;
		
		//convergence test
		if(max_point_translation_update<point_translation_small)
		{
			hasConverged=true;
			break;
		}
		
	}//end of for iter
}*/

void BundleAdjuster::BundleAdjustRobust2(int nb_iter)
{
	
	int nb_keyFrames=mCams.size();
	int nb_points=mPoints.size();
	
	float mLMLambda = 0.0001;//before LevMarConstant
	float mdLambdaFactor = 2.0;

	for(int iter=0;iter<nb_iter;iter++)
	{
		//std::cout<<"iteration : "<<iter <<std::endl;
		startIterRobust:
		
		//get Tukey factor
		std::vector<float> vdErrorSquared[nb_keyFrames];
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			if(mapPointsCam[2]>0)
			{
				//compute error with observed points (in meter in z=1 plane)
				Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
				Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
				Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
				vdErrorSquared[c].push_back(error.transpose()*error);
			}
		}
		float sigma_tukey[nb_keyFrames];
		for(int c=0;c<nb_keyFrames;c++)
		{
			sigma_tukey[c]=getSigmaSquared(vdErrorSquared[c]);
			if(sigma_tukey[c] < MINSIGMATUKEY)
				sigma_tukey[c] = MINSIGMATUKEY;
		}

		
		
		//points will be considered for optimisation only if they have two valid projection
		//that was kind of done already but Tukey function was not known yet
		//=> update taking robust function into account
		//std::cout<<"nb measure before outlier check : "<<mMeasures.size() <<std::endl;
		std::cout<<"nbCamsToUpdate : "<<nbCamsToUpdate <<std::endl;
		//std::cout<<"sigma_tukey : "<<sigma_tukey <<std::endl;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			  
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey[c]);
			//std::cout<<"kf :"<<c <<"  error = "<<(myCamera->m2PixProjJac()*(x_d-x_c)).transpose()<<std::endl;
			//std::cout<<"norm_reproj_error :"<<norm_reproj_error <<"  TukeyCoef = "<<TukeyCoef<<std::endl;
			if(mapPointsCam[2]<0 || mPoints[p].confidence*TukeyCoef==0)//measure is considered as outlier
			{
				//std::cout<<"reject measure from kf :"<<c <<std::endl;
				rejectMeasure(m);
				m--;
			}
		}
		
		//std::cout<<"nb measure after outlier check : "<<mMeasures.size() <<std::endl;
		std::cout<<"nbCamsToUpdate after: "<<nbCamsToUpdate <<std::endl;
		
		//check if BA is still possible
		if(nbCamsToUpdate==0)
		{
			OptimisePointPosition(nb_iter-iter);
			break;//has done point optim => goto end
		}
		
		//compute all non null elements of jacobien and reprojection error
		std::vector<kf_vis_jacobian> fJacobian;
		float residue_robust=0;
		float valid_wmeas=0;
		
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].pose;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].position;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey[c]);
			
			//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mapPointsCam[2]>0 && mPoints[p].confidence*TukeyCoef!=0 && (mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1))
			{
				kf_vis_jacobian jac;
				//get index of point to optim in list of point to optim
				jac.proj_error=x_d-x_c;
				jac.weight=TukeyCoef*mPoints[p].confidence;
				jac.pt_index=mPoints[p].id_optim;
				//get jacobian of error with respect to point position
				if(mPoints[p].id_optim!=-1)
					jac.de_dx=-myCamera->ProjectZ1_Jac_X(mapPointsCam)*pose.get_rotation();
				
				//get jacobien  of error with respect to variation of camera pose
				jac.cam_index=mCams[c].id_optim;
				if(mCams[c].id_optim!=-1)
					jac.de_dp=-myCamera->ProjectZ1_Jac_Dp(mapPointsCam);
				fJacobian.push_back(jac);
				
				//update residue to do levendberg marquardt check
				residue_robust+=mPoints[p].confidence*TukeyCoef*sqrt(norm_reproj_error);
				valid_wmeas+=mPoints[p].confidence*TukeyCoef;
				
				//if(isnan(jac.proj_error[0]))std::cout<<"nanery in error "<<std::endl;
			}
		}
		std::cout<<"residue_robust : "<<residue_robust <<std::endl;
		//std::cout<<"residue_robust norm: "<<residue_robust/valid_wmeas <<std::endl;
		

		//evaluate update using Gauss Newton
		//optim over points position + camera nb_keyFrames-1 pose = 3*nbPointsToUpdate + 6*(nb_keyFrames-1) dim
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
		
		
		std::vector<int> pt_gone_wrong;
		for(int i=0;i<nbPointsToUpdate;i++)
		{
			for(int j=0;j<3;j++)Hxx[i](j,j)=(1.+mLMLambda)*Hxx[i](j,j);
			//Eigen::FullPivLU<MatrixXf> luHxxi(Hxx[i]);
			//HxxInv[i]=luHxxi.inverse();
			HxxInv[i]=invMatrix3f(Hxx[i]);

			
			if(isnan(HxxInv[i](0,0)))
			{
				//find corresponding point and put in list not to be optimised anymore (if it is point at infinity then could still be used to locate camera pose)
				int id_pt=IdOptimToIdBA(i);
				pt_gone_wrong.push_back(id_pt);
				
			}
		}
		if(pt_gone_wrong.size()!=0)
		{
			coutErrMapping<<"\tBAR: strange, point hessian inversion gone wrong in iter "<< iter <<" = "<<pt_gone_wrong.size()<<endlErrMapping;
			for(int i=0;i<pt_gone_wrong.size();i++)
				doNotOptimisePoint(pt_gone_wrong[i]);
			
			//go back to start iteration loop
			goto startIterRobust;
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

		//update camera poses in temporary struct New

		float max_cam_translation_update=0;
		for(int c=0;c<nb_keyFrames;c++)
		{
			short id_opt=mCams[c].id_optim;
			if(id_opt!=-1)
			{
				
				VectorXf Dci=gain*Dc.segment(6*id_opt,6);
				mCams[c].poseNew=HomogeneousMatrix(Dci)*mCams[c].pose;
				
				if(isnan(Dci[0]))
				{
					//coutErrMapping<<"\tBA: strange, cam update"<< iter <<endlErrMapping;
					mCams[c].poseNew=mCams[c].pose;
					if(!goneWrong)
						std::cout<<"Dci["<<id_opt<<"] is nan" <<std::endl;
					goneWrong=true;
				}
				
				float translation_update=(Dci.segment(0,3)).squaredNorm();
				if(max_cam_translation_update<translation_update)
					max_cam_translation_update=translation_update;
			}
		}
		//std::cout<<"max_translation_update = "<<max_translation_update<<std::endl;


		//update point positions
		VectorXf Dx_all(3*nbPointsToUpdate);
		Dx_all=HxxInvJtp-HxxInvHxp*Dc;

		float max_point_translation_update=0;
		for(int i=0;i<nb_points;i++)
		{
			short id_opt=mPoints[i].id_optim;
			if(id_opt!=-1)
			{
				Vector3f Dx;for(int k=0;k<3;k++)Dx[k]=gain*Dx_all[3*id_opt+k];
				mPoints[i].positionNew=mPoints[i].position+Dx;
				
				if(isnan(Dx[0]))
				{
					//coutErrMapping<<"\tBA: strange, point update"<< iter <<endlErrMapping;
					mPoints[i].positionNew=mPoints[i].position;
					if(!goneWrong)
						std::cout<<"Dx["<<id_opt<<"] is nan" <<std::endl;
					goneWrong=true;
				}
				
				float translation_update=Dx.squaredNorm();
				if(max_point_translation_update<translation_update)
					max_point_translation_update=translation_update;
			}
		}
		
		//convergence test
		if(max_point_translation_update<point_translation_small && max_cam_translation_update<cam_translation_small)
		{
			hasConverged=true;
			break;
		}		
				
		//check if step was good
		float residue_robust_after=0;
		float valid_meas_after=0;
		int nb_outliers_new=0;
		for(int m=0;m<mMeasures.size();m++)
		{
			short &c=mMeasures[m].kf_id;
			short &p=mMeasures[m].pt_id;
			HomogeneousMatrix &pose=mCams[c].poseNew;
			
			//get 3D pose in camera frame
			Vector3f mapPointsCam=pose*mPoints[p].positionNew;	
			
			//compute error with observed points (in meter in z=1 plane)
			Vector2f x_d=mMeasures[m].coord;//desired projection= measurement
			Vector2f x_c=myCamera->ProjectZ1(mapPointsCam);//current projection				
			//Vector2f error=invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c);//error in pixels //WARNING has been modif, used to be in meters
			float norm_reproj_error=((invScaleLevel(mMeasures[m].lvl_meas)*myCamera->m2PixProjJac()*(x_d-x_c)).squaredNorm());
			float TukeyCoef=squareRootTukey(norm_reproj_error,sigma_tukey[c]);
			
			//TukeyCoef==0 should not happen anymore since outliers have been removed from list of measure just beforehand
			//compute jacobian only if point or cam linked to measure is to be optimised
			if(mPoints[p].confidence*TukeyCoef!=0 && (mPoints[p].id_optim!=-1 || mCams[c].id_optim!=1))
			{
				//residue to do levendberg marquardt check
				residue_robust_after+=mPoints[p].confidence*TukeyCoef*sqrt(norm_reproj_error);
				valid_meas_after+=mPoints[p].confidence*TukeyCoef;			}
			else
				nb_outliers_new++;
		}	
		
		coutMapping<<"residue before : "<<residue_robust<<"/"<<valid_wmeas<<" = "<<residue_robust/valid_wmeas<<"\t residue after : "<<residue_robust_after<<"/"<<valid_meas_after<<" = "<<residue_robust_after/valid_meas_after<<endlMapping;
					
		
		//if step was good confirm change and change lanbda
		if(residue_robust_after/valid_meas_after<residue_robust/valid_wmeas)
		{
			mdLambdaFactor = 2.0;
			mLMLambda *= 0.3;
			
			
			for(int c=0;c<nb_keyFrames;c++)
				mCams[c].pose=mCams[c].poseNew;
			for(int i=0;i<nb_points;i++)
				mPoints[i].position=mPoints[i].positionNew;
			
			if(nb_outliers_new!=0 && iter==nb_iter-1)
				//std::cout<<"There was "<<nb_outliers_new<<"new outliers"<<std::endl;
				coutColMapping<<"There was "<<nb_outliers_new<<"new outliers"<<endlColMapping;
			
			
		}
		//if not do not confirm change but change lanbda
		else
		{
			mLMLambda = mLMLambda * mdLambdaFactor;
			mdLambdaFactor = mdLambdaFactor * 2;
		}
#ifdef USE_OMP_C			
		if(canBeInterrupted)
		{
			omp_set_lock(lock_check_more_prior);
			if(*moreImportantStuffWaiting)
			{
			  omp_unset_lock(lock_check_more_prior);	
			  break; 
			}
			omp_unset_lock(lock_check_more_prior);	
		}
#endif	

		
	}//end of for iter



	//update all relative poses (dp not need to do it for each iteration)
	//std::vector<int> innerWin;for(int i=0;i<KeyFrameList.size();i++)innerWin.push_back(i);
	//updateRelativePoseInWin(innerWin);
	
}
