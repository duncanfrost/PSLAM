#include "PoseGraphOpt.h"

#define MIN_NB_MEASURE_PER_OPTIMISED_CAM 5

PoseGraphOptimiser::PoseGraphOptimiser(obMap *_myMap)
{
	myMap=_myMap;
	nbCamsToUpdate=0;
	hasConverged=false;
	cam_translation_small=1e-5;
	gain=DEFAULT_GAIN_GN;
}

void PoseGraphOptimiser::optimiseInnerWindow(std::vector<int> &innerWin,int nb_iter)
{
	std::cout<<"PoseGraphOptimiser optimiseInnerWindow"<<std::endl;
	std::vector<int> optimKF;
	optimKF=innerWin;
	//get outer fix windows:
	std::vector<int> fixedKF=myMap->getDirectNeigbors(optimKF);

	if(optimKF.size()!=0 && fixedKF.size()==0)
	{
		fixedKF.push_back(*optimKF.begin());
		optimKF.erase(optimKF.begin());
	}	
	if(optimKF.size()!=0)
	{
	
	//pass all the camera positions
	//std::cout<<"add KF"<<std::endl;
	for(int i=0;i<fixedKF.size();i++)
		addCamera(myMap->getKF(fixedKF[i])->getPose(),fixedKF[i],true);
	for(int i=0;i<optimKF.size();i++)
		addCamera(myMap->getKF(optimKF[i])->getPose(),optimKF[i],false);
	
	std::vector<int> allKF=optimKF;
	for(int i=0;i<fixedKF.size();i++)allKF.push_back(fixedKF[i]);

	
	for(int i=0;i<allKF.size();i++)
	{
		KeyFrame *ptKF=myMap->getKF(allKF[i]);
		for(int j=0;j<ptKF->getNbNeigbours();j++)
		{
			int id_neigbour=ptKF->getPtNeigbour(j)->neighboring_kf;
			//if this neigbot is in window then add edge
			if(std::find(allKF.begin(), allKF.end(), id_neigbour)!=allKF.end())
				if(allKF[i]<id_neigbour)//neigboring is defined in two direction, => if avoid duplicating edges
					//addEdge(allKF[i],id_neigbour,ptKF->getPtNeigbour(j)->relative_poses);
					//addEdge(id_neigbour,allKF[i],ptKF->getPtNeigbour(j)->relative_poses);
					//addEdge(id_neigbour,allKF[i],ptKF->getPtNeigbour(j)->relative_poses,ptKF->getPtNeigbour(j)->InformationMatrixPose);
					addEdge(id_neigbour,allKF[i],ptKF->getPtNeigbour(j)->relative_scale,ptKF->getPtNeigbour(j)->Informationscale,ptKF->getPtNeigbour(j)->relative_poses,ptKF->getPtNeigbour(j)->InformationMatrixPose);

		}
	}
	//optimise(nb_iter);
	optimiseSim3(nb_iter);
	
	//update poses
	for(int i=0;i<optimKF.size();i++)
		myMap->getKF(optimKF[i])->setPose(getUpdatedCamPose(optimKF[i]));
	
	//update scale if done on Sim3
	for(int k=0;k<optimKF.size();k++)
	{
		KeyFrame &KFc=*myMap->getKF(optimKF[k]);
		float scale_change=getUpdatedCamScale(optimKF[k]);
		std::cout<<"scale_change kf["<<optimKF[k]<<"]= "<<scale_change<<std::endl; 
		
		//rescale features
		for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
			KFc.getPtLocalBestFeatures(f)->depthInRef*=scale_change;
		
		//rescale best stereo
	  	HomogeneousMatrix relBestPose= KFc.getBestRelPose();
		relBestPose.set_translation(relBestPose.get_translation()*scale_change);
		KFc.setBestRelPose(relBestPose);
	}
	//update relative pose constraints
	for(int i=0;i<allKF.size();i++)
	{
		KeyFrame &KFc=*myMap->getKF(allKF[i]);
		for(int j=0;j<KFc.getNbNeigbours();j++)
		{
			int id_neigbour=KFc.getPtNeigbour(j)->neighboring_kf;
			KeyFrame &KFn=*myMap->getKF(id_neigbour);
			if(std::find(allKF.begin(), allKF.end(), id_neigbour)!=allKF.end())
			{
				NeigbourKFNew &neigbor=*KFc.getPtNeigbour(j);
				float scale_changec=getUpdatedCamScale(allKF[i]);
				float scale_changen=getUpdatedCamScale(neigbor.neighboring_kf);
				neigbor.relative_scale=scale_changec*neigbor.relative_scale/scale_changen;
				neigbor.relative_poses.set_translation(scale_changec*neigbor.relative_poses.get_translation());
			}
		}
	}
		
	
	//update mapPoint positions
	for(int k=0;k<optimKF.size();k++)
	{
		//get through all features and  if matched update corresponding point position
		KeyFrame &KFc=*myMap->getKF(optimKF[k]);
		for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
		{
			uptoscaleFeature &feat=*KFc.getPtLocalBestFeatures(f);
			if(feat.matched)
			{
				MapPoint &point=*feat.ptKForigin->getPtMapPoint(feat.idPoint);
				Vector3f newPos(0,0,0);
				float weightCurrent=0;
				for(int v=0;v<point.nbViews();v++)
				{
					KeyFrame &KFv=*myMap->getKF(point.getView(v));
					int id_feat_v=KFv.indexCandidateFeatureFromVisoId(point.getI1p(v));
					uptoscaleFeature &feat_v=*KFv.getPtLocalBestFeatures(id_feat_v);
					
					Vector3f localCoordn=toHomogeneous(feat_v.posRef)*feat_v.depthInRef;
					Vector3f AbsoluteCoordn=KFv.getPose().inverse()*localCoordn;
					//do average of coordinates of views weighted by scoreLocalFeature
					float weight_feat=feat_v.scoreFundamentalOrigin;
					weightCurrent+=weight_feat;
					newPos+=weight_feat*AbsoluteCoordn;
				}
				point.updatePosition(newPos/weightCurrent);
			}
		}
	}
	
	}
}

void PoseGraphOptimiser::addCamera(HomogeneousMatrix _pose,int _id_main,bool _fixed)
{
	//coutMapping<<"Add Camera "<<_id_main<<endlMapping;
	PGOcamera newCam;
	newCam.pose=_pose;
	newCam.id_main=_id_main;
	newCam.scale=1.;
	if(_fixed)
		newCam.id_optim=-1;
	else
	{
		newCam.id_optim=nbCamsToUpdate;
		nbCamsToUpdate++;
	}
	mCams.push_back(newCam);
}

void PoseGraphOptimiser::addEdge(int _id1,int _id2,HomogeneousMatrix pose12)
{
	PGOedge newEdge;
	newEdge.id_cam1=-1;
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id1)
		{
			newEdge.id_cam1=c;
			break;
		}
		
	newEdge.id_cam2=-1;
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id2)
		{
			newEdge.id_cam2=c;
			break;
		}
	newEdge.pose12=pose12;
	
	if(newEdge.id_cam1!=-1 && newEdge.id_cam2!=-1)
		mEdges.push_back(newEdge);
}
void PoseGraphOptimiser::addEdge(int _id1,int _id2,HomogeneousMatrix pose12,MatrixXf _InfoMatrix)
{
	PGOedge newEdge;
	newEdge.id_cam1=-1;
	
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id1)
		{
			newEdge.id_cam1=c;
			break;
		}
		
	newEdge.id_cam2=-1;
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id2)
		{
			newEdge.id_cam2=c;
			break;
		}
	
	if(newEdge.id_cam1!=-1 && newEdge.id_cam2!=-1)
	{
		newEdge.pose12=pose12;
		newEdge.InfoMatrix=_InfoMatrix;
		mEdges.push_back(newEdge);
	}
}

void PoseGraphOptimiser::addEdge(int _id1,int _id2,float scale12,float _infoScale,HomogeneousMatrix pose12,MatrixXf _InfoMatrix)
{
	PGOedge newEdge;
	newEdge.id_cam1=-1;
	
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id1)
		{
			newEdge.id_cam1=c;
			break;
		}
		
	newEdge.id_cam2=-1;
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==_id2)
		{
			newEdge.id_cam2=c;
			break;
		}
	
	if(newEdge.id_cam1!=-1 && newEdge.id_cam2!=-1)
	{
		newEdge.pose12=pose12;
		newEdge.InfoMatrix=_InfoMatrix;
		newEdge.scale12=scale12;
		newEdge.InfoScale=_infoScale;
		mEdges.push_back(newEdge);
	}
}



void PoseGraphOptimiser::optimise(int nb_iter)
{
	VectorXf Jacobian(6*nbCamsToUpdate);
	MatrixXf Hessian(6*nbCamsToUpdate,6*nbCamsToUpdate);
	if(nbCamsToUpdate<mCams.size())//need at least one fixed camera
	for(int iter=0;iter<nb_iter;iter++)
	{
		Jacobian.setZero();
		Hessian.setZero();
		for(int i=0;i<mEdges.size();i++)
		{
			int id_first_kf=mEdges[i].id_cam1;
			int id_second_kf=mEdges[i].id_cam2;
			
			//error to minimize = LogSe3(T01*T0*T1.inv)
			HomogeneousMatrix &relPose=mEdges[i].pose12;
			HomogeneousMatrix HError=relPose*mCams[id_first_kf].pose*mCams[id_second_kf].pose.inverse();
			VectorXf logError=HError.get_p();
			Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logError[j];
			Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logError[j+3];

			MatrixXf M1(6,6);
			M1.block(0,0,3,3)=-GetSkew(Dw);
			M1.block(3,3,3,3)=-GetSkew(Dw);
			M1.block(0,3,3,3)=-GetSkew(Dt);
			M1.block(3,0,3,3).setZero();
			
			MatrixXf M2(6,6);
			M2.block(0,0,3,3)=relPose.get_rotation();
			M2.block(3,3,3,3)=relPose.get_rotation();
			M2.block(0,3,3,3)=GetSkew(relPose.get_translation())*relPose.get_rotation();
			M2.block(3,0,3,3).setZero();
			
			//derivative wrt T0:
			MatrixXf dHErr_dp0=(MatrixXf::Identity(6,6)+0.5*M1)*M2;
			
			//derivative wrt T1:
			MatrixXf dHErr_dp1=-(MatrixXf::Identity(6,6)-0.5*M1);
			
			int id_opt_first_kf=mCams[mEdges[i].id_cam1].id_optim;
			int id_opt_second_kf=mCams[mEdges[i].id_cam2].id_optim;
			
			
			if(id_opt_first_kf!=-1)//set 0 as fix
			{
				Jacobian.segment(6*id_opt_first_kf,6)+=dHErr_dp0.transpose()*logError;
				Hessian.block(6*id_opt_first_kf,6*id_opt_first_kf,6,6)+=dHErr_dp0.transpose()*dHErr_dp0;
			}
			
			if(id_opt_second_kf!=-1)
			{
				Jacobian.segment(6*id_opt_second_kf,6)+=dHErr_dp1.transpose()*logError;
				Hessian.block(6*id_opt_second_kf,6*id_opt_second_kf,6,6)+=dHErr_dp1.transpose()*dHErr_dp1;			
			}

			if(id_opt_first_kf!=-1 && id_opt_second_kf!=-1)
			{
				Hessian.block(6*id_opt_second_kf,6*id_opt_first_kf,6,6)+=dHErr_dp1.transpose()*dHErr_dp0;
				Hessian.block(6*id_opt_first_kf,6*id_opt_second_kf,6,6)+=dHErr_dp0.transpose()*dHErr_dp1;
			}
		}

		

		Hessian.ldlt().solveInPlace(Jacobian);
		VectorXf Dp=-gain*Jacobian;
		
		
		float max_cam_translation_update=0;	
		for(int i=0;i<mCams.size();i++)
		{
			if(mCams[i].id_optim!=-1)
			{
				VectorXf tDp=Dp.segment(6*mCams[i].id_optim,6);
				mCams[i].pose=HomogeneousMatrix(tDp)*mCams[i].pose;
				
				float translation_update=(tDp.segment(0,3)).squaredNorm();
				if(max_cam_translation_update<translation_update)
					max_cam_translation_update=translation_update;
			}
		}
		//std::cout<<"max_cam_translation_update = "<<max_cam_translation_update<<std::endl;
		if(max_cam_translation_update<cam_translation_small)
		{
			hasConverged=true;
			break;
		}
	}
}

void PoseGraphOptimiser::optimiseSim3(int nb_iter)
{
	//will do optim changing temporary scale value, local reconstruction will be scales at the end 
	
	VectorXf Jacobian(7*nbCamsToUpdate);
	MatrixXf Hessian(7*nbCamsToUpdate,7*nbCamsToUpdate);
	if(nbCamsToUpdate<mCams.size())//need at least one fixed camera
	for(int iter=0;iter<nb_iter;iter++)
	{
		Jacobian.setZero();
		Hessian.setZero();
		for(int i=0;i<mEdges.size();i++)
		{
			int id_first_kf=mEdges[i].id_cam1;
			int id_second_kf=mEdges[i].id_cam2;
			
			//error to minimize = LogSe3(T01*T0*T1.inv)
			float scaleError=log(mEdges[i].scale12*mCams[id_second_kf].scale/mCams[id_first_kf].scale);
			//std::cout<<"mCams["<<i<<"].scale = "<<mCams[i].scale<<std::endl;
			//std::cout<<"mCams["<<i<<"].scale = "<<mCams[i].scale<<std::endl;
			
			HomogeneousMatrix relPose=mEdges[i].pose12;
			//goal is to have it equal to mCams[id_first_kf].pose*mCams[id_second_kf].pose.inverse()
			//relative translation constraint has been computed using init scale=1 of each camera, if scale of cam changes then this translation is scaled with it
			relPose.set_translation(mCams[id_first_kf].scale*relPose.get_translation());
			
			HomogeneousMatrix HError=relPose*mCams[id_first_kf].pose*mCams[id_second_kf].pose.inverse();
			VectorXf logError=HError.get_p();
			Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logError[j];
			Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logError[j+3];

			MatrixXf M1(6,6);
			M1.block(0,0,3,3)=-GetSkew(Dw);
			M1.block(3,3,3,3)=-GetSkew(Dw);
			M1.block(0,3,3,3)=-GetSkew(Dt);
			M1.block(3,0,3,3).setZero();
			
			MatrixXf M2(6,6);
			M2.block(0,0,3,3)=relPose.get_rotation();
			M2.block(3,3,3,3)=relPose.get_rotation();
			M2.block(0,3,3,3)=GetSkew(relPose.get_translation())*relPose.get_rotation();
			M2.block(3,0,3,3).setZero();
			
			//derivative wrt T0:
			MatrixXf dHErr_dp0=(MatrixXf::Identity(6,6)+0.5*M1)*M2;
			
			//derivative wrt T1:
			MatrixXf dHErr_dp1=-(MatrixXf::Identity(6,6)-0.5*M1);
			
			int id_opt_first_kf=mCams[mEdges[i].id_cam1].id_optim;
			int id_opt_second_kf=mCams[mEdges[i].id_cam2].id_optim;
			
			
			if(id_opt_first_kf!=-1)//set 0 as fix
			{
				Jacobian[7*id_opt_first_kf]+=-scaleError/mCams[id_first_kf].scale;
				Hessian(7*id_opt_first_kf,7*id_opt_first_kf)+=1./(mCams[id_first_kf].scale*mCams[id_first_kf].scale);
				Jacobian.segment(7*id_opt_first_kf+1,6)+=dHErr_dp0.transpose()*logError;
				Hessian.block(7*id_opt_first_kf+1,7*id_opt_first_kf+1,6,6)+=dHErr_dp0.transpose()*dHErr_dp0;
			}
			
			if(id_opt_second_kf!=-1)
			{
				Jacobian[7*id_opt_second_kf]+=scaleError/mCams[id_second_kf].scale;
				Hessian(7*id_opt_second_kf,7*id_opt_second_kf)+=1./(mCams[id_second_kf].scale*mCams[id_second_kf].scale);
				Jacobian.segment(7*id_opt_second_kf+1,6)+=dHErr_dp1.transpose()*logError;
				Hessian.block(7*id_opt_second_kf+1,7*id_opt_second_kf+1,6,6)+=dHErr_dp1.transpose()*dHErr_dp1;			
			}

			if(id_opt_first_kf!=-1 && id_opt_second_kf!=-1)
			{
				Hessian(7*id_opt_second_kf,7*id_opt_first_kf)+=-1./(mCams[id_first_kf].scale*mCams[id_second_kf].scale);
				Hessian(7*id_opt_first_kf,7*id_opt_second_kf)+=-1./(mCams[id_first_kf].scale*mCams[id_second_kf].scale);
				Hessian.block(7*id_opt_second_kf+1,7*id_opt_first_kf+1,6,6)+=dHErr_dp1.transpose()*dHErr_dp0;
				Hessian.block(7*id_opt_first_kf+1,7*id_opt_second_kf+1,6,6)+=dHErr_dp0.transpose()*dHErr_dp1;
			}
		}

		

		Hessian.ldlt().solveInPlace(Jacobian);
		VectorXf Dp=-gain*Jacobian;
		
		
		float max_cam_translation_update=0;	
		for(int i=0;i<mCams.size();i++)
		{
			if(mCams[i].id_optim!=-1)
			{
				float tDs=Dp[7*mCams[i].id_optim];
				mCams[i].scale=mCams[i].scale+tDs;
				//std::cout<<"mCams["<<i<<"].scale = "<<mCams[i].scale<<std::endl;
				VectorXf tDp=Dp.segment(7*mCams[i].id_optim+1,6);
				mCams[i].pose=HomogeneousMatrix(tDp)*mCams[i].pose;
				
				float translation_update=(tDp.segment(0,3)).squaredNorm();
				if(max_cam_translation_update<translation_update)
					max_cam_translation_update=translation_update;
			}
		}
		//std::cout<<"max_cam_translation_update = "<<max_cam_translation_update<<std::endl;
		/*if(max_cam_translation_update<cam_translation_small)
		{
			hasConverged=true;
			break;
		}*/
	}
}


HomogeneousMatrix PoseGraphOptimiser::getUpdatedCamPose(int id_main)
{
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==id_main)return mCams[c].pose;
		
	std::cerr<<"BundleAjuster::getUpdatedCamPose> cam not found"<<endlMapping;
	return HomogeneousMatrix();//should not happen
}

float PoseGraphOptimiser::getUpdatedCamScale(int id_main)
{
	for(int c=0;c<mCams.size();c++)
		if(mCams[c].id_main==id_main)return mCams[c].scale;
		
	//std::cerr<<"BundleAjuster::getUpdatedCamScale> cam not found"<<endlMapping;
	return 1;//should not happen
}

