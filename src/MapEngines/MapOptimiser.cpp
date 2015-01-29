#include "MapOptimiser.h"

#define MIN_NB_MEASURE_PER_OPTIMISED_CAM 5

MapOptimiser::MapOptimiser(obMap *_myMap)
{
	myMap=_myMap;
	gain=DEFAULT_GAIN_GN;
}


struct kf_loc_jacobian
{
	Vector3f pos_error;//projection error
	int opt_id;//index of kf in optimisation list
	MatrixXf de_dp;//jacobian wrt kf position
	Vector3f de_ds;//jacobian wrt kf scale
	float weight;
};
struct kf_loc_jacobian2
{
	Vector3f pos_error;//projection error
	int opt_idc;//index of kf in optimisation list
	MatrixXf de_dpc;//jacobian wrt kf position
	Vector3f de_dsc;//jacobian wrt kf scale
	int opt_idn;//index of kf in optimisation list
	MatrixXf de_dpn;//jacobian wrt kf position
	Vector3f de_dsn;//jacobian wrt kf scale
	float weight;
};
void MapOptimiser::optimiseInnerWindowRobust(std::vector<int> &_innerWindowKFs,int nb_iter)
{
	//coutGreen<<"###############################################"<<endlGreen;
	coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  "<<std::endl;
	//need one KF to be fixed
	//if empty then need to fix one frame, take first one out
	if(outerWindow.size()==0)
	{
		outerWindow.push_back(*innerWindowKFs.begin());
		innerWindowKFs.erase(innerWindowKFs.begin());
	}

  	if(_innerWindowKFs.size()!=0)
	{
	//union of inner and outer
	std::vector<int> FullWindow=innerWindowKFs;
	for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
	
	//check windows:
	//std::cout<<"FullWindow = "<<std::endl;
	//for(int i=0;i<FullWindow.size();i++)std::cout<<FullWindow[i]<<"  "<<std::endl;
	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	//std::cout<<std::endl;
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	//std::cout<<std::endl<<std::endl;
	
	bool verb_BA=false;
	int nbOptimKf=innerWindowKFs.size();
	
	Camera *myCam=myMap->getCamera();
		
	if(nbOptimKf>0)
	{

	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		HomogeneousMatrix22 pose_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)pose_kfs[k]=myMap->getKF(innerWindowKFs[k])->getPose();
		
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		  
		
		if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
		for(int iter=0;iter<nb_iter;iter++)
		{
			//if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			//if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			std::cout<<"\t###############################################"<<std::endl;
			std::cout<<"\tIteration "<<iter<<std::endl;
			//need to put in minimisation all the points linked by features of innerWindow and outerWindow
			//that should be all the points defined in this keyframes
			std::vector<kf_loc_jacobian2> list_jacobian_scale;
			
			//residueEssential just to check if error goes down
			float residue=0;
			float nb_residue=0;
			for(int k=0;k<FullWindow.size();k++)
			{
				int &idKF=FullWindow[k];
				KeyFrame &KFc=*myMap->getKF(idKF);
				int id_opt_c=LUTidKFtoInnerWindow[idKF];
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb points : "<<KFc.getNbMapPoint()<<std::endl;
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb Feat : "<<KFc.getNbLocalBestFeatures()<<std::endl;
				
				//get neigbors
				for(int n=0;n<KFc.getNbNeigbours();n++)
				{
					NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
					KeyFrame &KFn=*myMap->getKF(neigbor.neighboring_kf);
					int id_opt_n=LUTidKFtoInnerWindow[neigbor.neighboring_kf];
					if(neigbor.neighboring_kf<idKF)//just consider one way, other way should be same error
					//if(neigbor.neighboring_kf>idKF)//just consider one way, other way should be same error
					{
						//if(verb_BA)std::cout<<"\tuse edge "<<idKF<<" <=> "<<neigbor.neighboring_kf<<std::endl;
						//compute essential matrix error; ie x_i *E_ij *x_j

					  
					  
						//compute essential matrix
						HomogeneousMatrix22 pose_c=KFc.getPose();
						HomogeneousMatrix22 pose_n=KFn.getPose();
						
						//HomogeneousMatrix relPose=neigbor.relative_poses;
						HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
						
						//get TukeyScalar
						std::vector<float> vdErrorSquared;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							p_match &mMatch=neigbor.matches[m];
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
							//get current estimated scales
							float scale_c;
							if(id_opt_c==-1)scale_c=1;
							else scale_c=scales_kfs[id_opt_c];
								
							float scale_n;
							if(id_opt_n==-1)scale_n=1;
							else scale_n=scales_kfs[id_opt_n];
								
							uptoscaleFeature *feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
							uptoscaleFeature *feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
							//float add_scales=(scale_c+scale_n);
							//Vector3f scaleError=(1./add_scales)*(scale_c*feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates()));
							Vector3f scaleError=feat_c->getLocalCoordinates()-(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
							vdErrorSquared.push_back(scaleError.squaredNorm());
							}
						}
						
						float sigma_tukey=0;
						if(vdErrorSquared.size()!=0)
							sigma_tukey=getSigmaSquared(vdErrorSquared);
						
						if(verb_BA)std::cout<<"\tnb matches = "<<neigbor.matches.size()<<std::endl;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							p_match &mMatch=neigbor.matches[m];
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
							//get current estimated scales
							float scale_c;
							if(id_opt_c==-1)scale_c=1;
							else scale_c=scales_kfs[id_opt_c];
								
							float scale_n;
							if(id_opt_n==-1)scale_n=1;
							else scale_n=scales_kfs[id_opt_n];
								
							uptoscaleFeature *feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
							uptoscaleFeature *feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
							//float add_scales=(scale_c+scale_n);
							//Vector3f scaleError=(1./add_scales)*(scale_c*feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates()));
							Vector3f scaleError=feat_c->getLocalCoordinates()-(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
							float TukeyCoef=squareRootTukey(scaleError.squaredNorm(),sigma_tukey);
							
							residue+=scaleError.squaredNorm();
							nb_residue+=TukeyCoef;
							
							if(TukeyCoef>0)
							{
								//deriv with respect to relPose:
								MatrixXf jac_dp(3,6);
								jac_dp.block(0,0,3,3)=-Matrix3f::Identity();
								//jac_dp.block(0,3,3,3)=GetSkew(relPose*(scale_n*feat_n->getLocalCoordinates()));
								jac_dp.block(0,3,3,3)=(1./scale_c)*GetSkew(relPose*(scale_n*feat_n->getLocalCoordinates()));
								//jac_dp=jac_dp/add_scales;
								
								//to get from relPose to deriv wrt pose_c and pose_n
								MatrixXf M1(6,6);
								VectorXf logRelPose=relPose.get_p();
								Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logRelPose[j];
								Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logRelPose[j+3];
								M1.block(0,0,3,3)=-GetSkew(Dw);
								M1.block(3,3,3,3)=-GetSkew(Dw);
								M1.block(0,3,3,3)=-GetSkew(Dt);
								M1.block(3,0,3,3).setZero();
								
								kf_loc_jacobian2 mJacScale;
								mJacScale.pos_error=scaleError;
								mJacScale.weight=TukeyCoef;
								
								mJacScale.opt_idc=id_opt_c;
								MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6)+0.5*M1);
								mJacScale.de_dpc=jac_dp*dHErr_dpc;
								//mJacScale.de_dsc=feat_c->getLocalCoordinates();
								//mJacScale.de_dsc=scale_n/(add_scales*add_scales)*feat_c->getLocalCoordinates();
								mJacScale.de_dsc=1./(scale_c*scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
							
							
								mJacScale.opt_idn=id_opt_n;
								MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1);
								mJacScale.de_dpn=jac_dp*dHErr_dpn;
								//mJacScale.de_dsn=-relPose.get_rotation()*  feat_n->getLocalCoordinates();
								//mJacScale.de_dsn=-scale_c/(add_scales*add_scales)*relPose.get_rotation()*  feat_n->getLocalCoordinates()+relPose.get_translation()/(add_scales*add_scales);
								mJacScale.de_dsn=-(1./scale_c)*relPose.get_rotation()*  feat_n->getLocalCoordinates();
								list_jacobian_scale.push_back(mJacScale);
							}
							}

						}
					}
				}
								
			}
			
			std::cout<<"\nresidue = "<<residue/nb_residue<<std::endl;

			
			//accumulate jacobians
			int nb_params=7*nbOptimKf;//scale + translation + rotation
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			
			//std::cout<<"Update Matrices using jacobianScale"<<std::endl;
			for(int j=0;j<list_jacobian_scale.size();j++)
			{
				kf_loc_jacobian2 &fJacobian=list_jacobian_scale[j];
				if(fJacobian.opt_idc!=-1)
				{
				Jte[7*fJacobian.opt_idc]+=fJacobian.weight * fJacobian.de_dsc.transpose()*fJacobian.pos_error;
				Jte.segment(7*fJacobian.opt_idc+1,6)+=fJacobian.weight * fJacobian.de_dpc.transpose()*fJacobian.pos_error;
				H(7*fJacobian.opt_idc,7*fJacobian.opt_idc)+=fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsc;				
				H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpc;				
				
				}
				
				if(fJacobian.opt_idn!=-1)
				{
				Jte[7*fJacobian.opt_idn]+=fJacobian.weight * fJacobian.de_dsn.transpose()*fJacobian.pos_error;
				Jte.segment(7*fJacobian.opt_idn+1,6)+=fJacobian.weight * fJacobian.de_dpn.transpose()*fJacobian.pos_error;
				H(7*fJacobian.opt_idn,7*fJacobian.opt_idn)+=fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsn;				
				H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpn;				
				}
				
				
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1)
				{
				H(7*fJacobian.opt_idc,7*fJacobian.opt_idn)+=fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsn;				
				H(7*fJacobian.opt_idn,7*fJacobian.opt_idc)+=fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsc;				
				H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpc;				
				H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpn;
				}

			}

			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
				VectorXf Dp(nb_params);
				Dp=-gain*(invH*Jte);
				
				if(verb_BA)std::cout<<"\tapply update"<<std::endl;
				//apply update to kf
				for(int j=0;j<nbOptimKf;j++)
				{
					if(verb_BA)std::cout<<"\tapply update scale"<<std::endl;
					float ds=Dp[7*j];
					//scales_kfs[j]=(1.+ds)*scales_kfs[j];
					scales_kfs[j]=ds+scales_kfs[j];
					VectorXf dp(Dp.segment(7*j+1,6));
					if(verb_BA)std::cout<<dp<<std::endl;
					pose_kfs[j]=HomogeneousMatrix22(dp)* pose_kfs[j];
				}
				
				
				if(verb_BA)std::cout<<"\tupdate map with new estimates"<<std::endl;
				//update KFs with new scale and pose //make it faster to compute point position of do that here in loop
				for(int k=0;k<nbOptimKf;k++)
				{
					if(verb_BA)std::cout<<"\t\tupdate "<<k<<std::endl;
					int &idKF=innerWindowKFs[k];
					KeyFrame &KFc=*myMap->getKF(idKF);
					
					VectorXf dp(Dp.segment(7*k+1,6));
					if(verb_BA)std::cout<<"\t\tpose up "<<HomogeneousMatrix22(dp)<<std::endl;
					if(verb_BA)std::cout<<"\t\tscale up "<<scales_kfs[k]<<std::endl;
					KFc.setPose(pose_kfs[k]);
					for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
						KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
					
					HomogeneousMatrix22 relPose= KFc.getRelativePose();
					relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
					KFc.setRelativePose(relPose);
					
					HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
					relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
					KFc.setBestRelPose(relBestPose);
					
					scales_kfs[k]=1;
				}
				
				if(verb_BA)std::cout<<"\tupdate point positions"<<std::endl;
				//update position of points
				for(int k=0;k<FullWindow.size();k++)
				{
					int &idKF=FullWindow[k];
					//std::cout<<"\tidKF "<<idKF<<std::endl;
					KeyFrame &KFc=*myMap->getKF(idKF);
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
							KeyFrame &KFv=*myMap->getKF(kfView);
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
			
		}
	}
}
	
}

void MapOptimiser::optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter)
{
	//coutGreen<<"###############################################"<<endlGreen;
	coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  "<<std::endl;
	//need one KF to be fixed
	//if empty then need to fix one frame, take first one out
	if(outerWindow.size()==0)
	{
		outerWindow.push_back(*innerWindowKFs.begin());
		innerWindowKFs.erase(innerWindowKFs.begin());
	}

  	if(_innerWindowKFs.size()!=0)
	{
	//union of inner and outer
	std::vector<int> FullWindow=innerWindowKFs;
	for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
	
	//check windows:
	//std::cout<<"FullWindow = "<<std::endl;
	//for(int i=0;i<FullWindow.size();i++)std::cout<<FullWindow[i]<<"  "<<std::endl;
	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	//std::cout<<std::endl;
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	//std::cout<<std::endl<<std::endl;
	
	bool verb_BA=false;
	int nbOptimKf=innerWindowKFs.size();
	
	Camera *myCam=myMap->getCamera();
		
	if(nbOptimKf>0)
	{

	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		HomogeneousMatrix22 pose_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)pose_kfs[k]=myMap->getKF(innerWindowKFs[k])->getPose();
		
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		  
		
		if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
		for(int iter=0;iter<nb_iter;iter++)
		{
			//if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			//if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			std::cout<<"\t###############################################"<<std::endl;
			std::cout<<"\tIteration "<<iter<<std::endl;
			//need to put in minimisation all the points linked by features of innerWindow and outerWindow
			//that should be all the points defined in this keyframes
			std::vector<kf_loc_jacobian2> list_jacobian_scale;
			
			float residue=0;
			float nb_residue=0;
			//residueEssential just to check if error goes down
			for(int k=0;k<FullWindow.size();k++)
			{
				int &idKF=FullWindow[k];
				KeyFrame &KFc=*myMap->getKF(idKF);
				int id_opt_c=LUTidKFtoInnerWindow[idKF];
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb points : "<<KFc.getNbMapPoint()<<std::endl;
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb Feat : "<<KFc.getNbLocalBestFeatures()<<std::endl;
				
				//get neigbors
				for(int n=0;n<KFc.getNbNeigbours();n++)
				{
					NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
					KeyFrame &KFn=*myMap->getKF(neigbor.neighboring_kf);
					int id_opt_n=LUTidKFtoInnerWindow[neigbor.neighboring_kf];
					if(neigbor.neighboring_kf<idKF)//just consider one way, other way should be same error
					//if(neigbor.neighboring_kf>idKF)//just consider one way, other way should be same error
					{
						//if(verb_BA)std::cout<<"\tuse edge "<<idKF<<" <=> "<<neigbor.neighboring_kf<<std::endl;
						//compute essential matrix error; ie x_i *E_ij *x_j

						//compute essential matrix
						HomogeneousMatrix22 pose_c=KFc.getPose();
						if(id_opt_c!=-1)pose_c=pose_kfs[id_opt_c];
						HomogeneousMatrix22 pose_n=KFn.getPose();
						if(id_opt_n!=-1)pose_n=pose_kfs[id_opt_n];
						
						//HomogeneousMatrix relPose=neigbor.relative_poses;
						HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
						
						if(verb_BA)std::cout<<"\tnb matches = "<<neigbor.matches.size()<<std::endl;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							p_match &mMatch=neigbor.matches[m];
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
							//get current estimated scales
							float scale_c;
							if(id_opt_c==-1)scale_c=1;
							else scale_c=scales_kfs[id_opt_c];
								
							float scale_n;
							if(id_opt_n==-1)scale_n=1;
							else scale_n=scales_kfs[id_opt_n];
								
							uptoscaleFeature *feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
							uptoscaleFeature *feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
							//float add_scales=(scale_c+scale_n);
							//Vector3f scaleError=scale_c*feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates());
							//Vector3f scaleError=(1./add_scales)*(scale_c*feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates()));
							Vector3f scaleError=feat_c->getLocalCoordinates()-(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
							residue+=scaleError.squaredNorm();
							nb_residue++;
							
							
							//deriv with respect to relPose:
							MatrixXf jac_dp(3,6);
							jac_dp.block(0,0,3,3)=-Matrix3f::Identity();
							//jac_dp.block(0,3,3,3)=GetSkew(relPose*(scale_n*feat_n->getLocalCoordinates()));
							jac_dp.block(0,3,3,3)=(1./scale_c)*GetSkew(relPose*(scale_n*feat_n->getLocalCoordinates()));
							//jac_dp=jac_dp/add_scales;
							
							//to get from relPose to deriv wrt pose_c and pose_n
							MatrixXf M1(6,6);
							VectorXf logRelPose=relPose.get_p();
							Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logRelPose[j];
							Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logRelPose[j+3];
							M1.block(0,0,3,3)=-GetSkew(Dw);
							M1.block(3,3,3,3)=-GetSkew(Dw);
							M1.block(0,3,3,3)=-GetSkew(Dt);
							M1.block(3,0,3,3).setZero();
							
							kf_loc_jacobian2 mJacScale;
							mJacScale.pos_error=scaleError;
							mJacScale.weight=1.;
							
							mJacScale.opt_idc=id_opt_c;
							MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6)+0.5*M1);
							mJacScale.de_dpc=jac_dp*dHErr_dpc;
							//mJacScale.de_dsc=feat_c->getLocalCoordinates();
							//mJacScale.de_dsc=scale_n/(add_scales*add_scales)*feat_c->getLocalCoordinates();
							mJacScale.de_dsc=1./(scale_c*scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
							
							
							mJacScale.opt_idn=id_opt_n;
							MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1);
							mJacScale.de_dpn=jac_dp*dHErr_dpn;
							//mJacScale.de_dsn=-relPose.get_rotation()*  feat_n->getLocalCoordinates();
							//mJacScale.de_dsn=-scale_c/(add_scales*add_scales)*relPose.get_rotation()*  feat_n->getLocalCoordinates()+relPose.get_translation()/(add_scales*add_scales);
							mJacScale.de_dsn=-(1./scale_c)*relPose.get_rotation()*  feat_n->getLocalCoordinates();
							list_jacobian_scale.push_back(mJacScale);
							}

						}
					}
				}
								
			}
			
			//accumulate jacobians
			int nb_params=7*nbOptimKf;//scale + translation + rotation
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			
			//std::cout<<"Update Matrices using jacobianScale"<<std::endl;
			for(int j=0;j<list_jacobian_scale.size();j++)
			{
				kf_loc_jacobian2 &fJacobian=list_jacobian_scale[j];
				if(fJacobian.opt_idc!=-1)
				{
				Jte[7*fJacobian.opt_idc]+=fJacobian.weight * fJacobian.de_dsc.transpose()*fJacobian.pos_error;
				Jte.segment(7*fJacobian.opt_idc+1,6)+=fJacobian.weight * fJacobian.de_dpc.transpose()*fJacobian.pos_error;
				H(7*fJacobian.opt_idc,7*fJacobian.opt_idc)+=fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsc;				
				H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpc;				
				
				}
				
				if(fJacobian.opt_idn!=-1)
				{
				Jte[7*fJacobian.opt_idn]+=fJacobian.weight * fJacobian.de_dsn.transpose()*fJacobian.pos_error;
				Jte.segment(7*fJacobian.opt_idn+1,6)+=fJacobian.weight * fJacobian.de_dpn.transpose()*fJacobian.pos_error;
				H(7*fJacobian.opt_idn,7*fJacobian.opt_idn)+=fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsn;				
				H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpn;				
				}
				
				
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1)
				{
				H(7*fJacobian.opt_idc,7*fJacobian.opt_idn)+=fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsn;				
				H(7*fJacobian.opt_idn,7*fJacobian.opt_idc)+=fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsc;				
				H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpc;				
				H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpn;
				}

			}

			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
				VectorXf Dp(nb_params);
				Dp=-gain*(invH*Jte);
				
				if(verb_BA)std::cout<<"\tapply update"<<std::endl;
				//apply update to kf
				for(int j=0;j<nbOptimKf;j++)
				{
					if(verb_BA)std::cout<<"\tapply update scale"<<std::endl;
					float ds=Dp[7*j];
					//scales_kfs[j]=(1.+ds)*scales_kfs[j];
					scales_kfs[j]=ds+scales_kfs[j];
					VectorXf dp(Dp.segment(7*j+1,6));
					if(verb_BA)std::cout<<dp<<std::endl;
					pose_kfs[j]=HomogeneousMatrix22(dp)* pose_kfs[j];
				}
				
			}
		}
			
		if(verb_BA)std::cout<<"\tupdate map with new estimates"<<std::endl;
		//update KFs with new scale and pose //make it faster to compute point position of do that here in loop
		for(int k=0;k<nbOptimKf;k++)
		{
			if(verb_BA)std::cout<<"\t\tupdate "<<k<<std::endl;
			int &idKF=innerWindowKFs[k];
			KeyFrame &KFc=*myMap->getKF(idKF);
			
			if(verb_BA)std::cout<<"\t\tscale up "<<scales_kfs[k]<<std::endl;
			KFc.setPose(pose_kfs[k]);
			for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
				KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
			
			HomogeneousMatrix22 relPose= KFc.getRelativePose();
			relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
			KFc.setRelativePose(relPose);
			
			HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
			relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
			KFc.setBestRelPose(relBestPose);
			
			scales_kfs[k]=1;
		}
		
		if(verb_BA)std::cout<<"\tupdate point positions"<<std::endl;
		//update position of points
		for(int k=0;k<FullWindow.size();k++)
		{
			int &idKF=FullWindow[k];
			//std::cout<<"\tidKF "<<idKF<<std::endl;
			KeyFrame &KFc=*myMap->getKF(idKF);
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
					KeyFrame &KFv=*myMap->getKF(kfView);
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
}
	
}

void MapOptimiser::optimiseInnerWindow2(std::vector<int> &_innerWindowKFs,int nb_iter)
{
  	if(_innerWindowKFs.size()!=0)
	{
	coutGreen<<"###############################################"<<endlGreen;
	coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  "<<std::endl;
	//need one KF to be fixed
	//if empty then need to fix one frame, take first one out
	if(outerWindow.size()==0)
	{
		outerWindow.push_back(*innerWindowKFs.begin());
		innerWindowKFs.erase(innerWindowKFs.begin());
	}

	//union of inner and outer
	std::vector<int> FullWindow=innerWindowKFs;
	for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
	
	//check windows:
	//std::cout<<"FullWindow = "<<std::endl;
	//for(int i=0;i<FullWindow.size();i++)std::cout<<FullWindow[i]<<"  "<<std::endl;
	std::cout<<"innerWindowKFs = "<<std::endl;
	for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	std::cout<<std::endl;
	std::cout<<"outerWindow = "<<std::endl;
	for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	std::cout<<std::endl<<std::endl;
	
	bool verb_BA=false;
	int nbOptimKf=innerWindowKFs.size();
	
		
	if(nbOptimKf>0)
	{

	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		HomogeneousMatrix22 pose_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)pose_kfs[k]=myMap->getKF(innerWindowKFs[k])->getPose();
		
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		  
		
		if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
		for(int iter=0;iter<nb_iter;iter++)
		{
			//if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			//if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			std::cout<<"\t###############################################"<<std::endl;
			std::cout<<"\tIteration "<<iter<<std::endl;
			//need to put in minimisation all the points linked by features of innerWindow and outerWindow
			//that should be all the points defined in this keyframes
			std::vector<kf_loc_jacobian> list_jacobian;
			
			//residue just to check if error goes down
			float residue=0;
			int p_disp=0;//only for debugging
			for(int k=0;k<FullWindow.size();k++)
			{
				int &idKF=FullWindow[k];
				KeyFrame &KFc=*myMap->getKF(idKF);
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb points : "<<KFc.getNbMapPoint()<<std::endl;
				for(int p=0;p<KFc.getNbMapPoint();p++)
				{
					
					if(p<p_disp)std::cout<<"\t\tpoint "<<p<<std::endl;
					//now we have the point and each of its views
					//each view has to be used to compute the sum to minimise
					MapPoint *point=KFc.getPtMapPoint(p);
					if(point->isUsed())//if point is good
					{
						if(p<p_disp)std::cout<<"\t\tpointInCurrentKF = "<<point->getPosition().transpose()<<std::endl;

						for(int v=0;v<point->nbViews();v++)
						{
							if(p<p_disp)std::cout<<"\t\t\tview = "<<v<<std::endl;
							int kfView=point->getView(v);
							if(p<p_disp)std::cout<<"\t\t\tkfview = "<<kfView<<std::endl;
							int i1pView=point->getI1p(v);
							//get corresponding local feature:
							int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
							if(p<p_disp)std::cout<<"\t\t\tid_feat_v = "<<id_feat_v<<std::endl;
							if(id_feat_v!=-1)
							{
								uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
								int id_opt=LUTidKFtoInnerWindow[kfView];
								if(id_opt!=-1)//if corresponding keyframe is to be optimised
								{
									//compute error position point versus position local feature scaled 
									Vector3f pointInCurrentKF=pose_kfs[id_opt]*point->getPosition();
									Vector3f pointFromLocalFeat=scales_kfs[id_opt]*feat_v.getLocalCoordinates();
									Vector3f positionError= pose_kfs[id_opt]*point->getPosition()-pointFromLocalFeat;
									if(p<p_disp)
									{
									//std::cout<<"\t\t\tpointFromLocalFeat = "<<pointFromLocalFeat.transpose()<<std::endl;
									std::cout<<"\t\t\tkfview : "<<kfView<<" Coord view : "<<feat_v.getLocalCoordinates().transpose()<<" weight = "<<feat_v.scoreFundamentalOrigin<<std::endl;
									}
									
									residue+=positionError.squaredNorm();
									
									//create jacobian
									kf_loc_jacobian fJacobian;
									fJacobian.pos_error=positionError;
									fJacobian.opt_id=id_opt;
									
									//compute jacobians:
									//jacobian wrt kf position: consider that search for update dp to compose with pose_kfs[id_opt] to minise equation
									//=> derive H(dp)*pointInCurrentKF, need d H(dp)/d dp
									//d H(dp) is computed at dp=0
									//where H(dp) X = R X + t = (I+Ohmega(r))X+t
									//with Ohmega is skew symetric matric =  [0    -r_z   r_y ]
									//					  [r_z  0     -r_x]
									//					  [-r_y r_x    0  ]
									//=> R X = [x       - r_z * y + r_y* z
									//=>     = [r_z * x + y       - r_x* z
									//=>     = [-r_y * x+ r_x * y + z
									//=> d R X / d r = [0   z  -y]  = -[X]x
									//                 [-z  0   x] 
									//                 [y  -x   0] 
									
									//has to be 3by6: diff of pos error 3 by 6 dof
									MatrixXf de_dp_t(3,6);
									de_dp_t(0,0)=1; de_dp_t(0,1)=0; de_dp_t(0,2)=0;    de_dp_t(0,3)=0;                    de_dp_t(0,4)=pointInCurrentKF[2];   de_dp_t(0,5)=-pointInCurrentKF[1];      
									de_dp_t(1,0)=0; de_dp_t(1,1)=1; de_dp_t(1,2)=0;    de_dp_t(1,3)=-pointInCurrentKF[2]; de_dp_t(1,4)=0;                     de_dp_t(1,5)= pointInCurrentKF[0];      
									de_dp_t(2,0)=0; de_dp_t(2,1)=0; de_dp_t(2,2)=1;    de_dp_t(2,3)=pointInCurrentKF[1];  de_dp_t(2,4)=-pointInCurrentKF[0];  de_dp_t(2,5)=0;      
									
									fJacobian.de_dp=de_dp_t;
									
									
									//jacobian wrt kf scale
									//Xlocal2=(1+ds)s*Xlocal=>
									fJacobian.de_ds=-scales_kfs[id_opt]*feat_v.getLocalCoordinates();
									
									
									
									//can use weight for this measure as reconstruction angle point or fundamental matrix score of origin
									fJacobian.weight=feat_v.scoreFundamentalOrigin;
									
									list_jacobian.push_back(fJacobian);
								}
							}
							
						}
					}
				}				
			}
			coutRed<<"\tresitue BA = "<<residue<<endlRed;
			if(verb_BA)std::cout<<"\tresitue BA normalized= "<<residue/list_jacobian.size()<<std::endl;
			std::cout<<"\tlist_jacobian size = "<<list_jacobian.size()<<std::endl;
			
			//compute update from jacobians
			int nb_params=7*nbOptimKf;
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			if(verb_BA)std::cout<<"\taccumulate in J and H"<<std::endl;
			for(int j=0;j<list_jacobian.size();j++)
			{
				kf_loc_jacobian &fJacobian=list_jacobian[j];
				
				//scale first
				Jte[7*fJacobian.opt_id]+=fJacobian.weight * fJacobian.de_ds.transpose()*fJacobian.pos_error;
				//then kf position
				Jte.segment(7*fJacobian.opt_id+1,6)+=fJacobian.weight * fJacobian.de_dp.transpose()*fJacobian.pos_error;
				
				//Hessian:
				H(7*fJacobian.opt_id,7*fJacobian.opt_id)+=fJacobian.weight * fJacobian.de_ds.transpose() * fJacobian.de_ds;				
				H.block(7*fJacobian.opt_id+1,7*fJacobian.opt_id+1,6,6)+=fJacobian.weight * fJacobian.de_dp.transpose() * fJacobian.de_dp;				
			}
			
			if(verb_BA)std::cout<<"\tcompute update"<<std::endl;
			//std::cout<<"\tJte = "<<std::endl;
			//std::cout<<Jte<<std::endl;
			//std::cout<<"\tH = "<<std::endl;
			//std::cout<<H<<std::endl;
			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
			  
				VectorXf Dp(nb_params);
				Dp=-gain*(invH*Jte);
				
				if(verb_BA)std::cout<<"\tapply update"<<std::endl;
				//apply update to kf
				for(int j=0;j<nbOptimKf;j++)
				{
					if(verb_BA)std::cout<<"\tapply update scale"<<std::endl;
					float ds=Dp[7*j];
					scales_kfs[j]=(1.+ds)*scales_kfs[j];
					
					if(verb_BA)std::cout<<"\tapply update pose"<<std::endl;
					VectorXf dp(Dp.segment(7*j+1,6));
					pose_kfs[j]=HomogeneousMatrix22(dp)* pose_kfs[j];
				}
				
				
				if(verb_BA)std::cout<<"\tupdate map with new estimates"<<std::endl;
				//update KFs with new scale and pose //make it faster to compute point position of do that here in loop
				for(int k=0;k<nbOptimKf;k++)
				{
					if(verb_BA)std::cout<<"\t\tupdate "<<k<<std::endl;
					int &idKF=innerWindowKFs[k];
					KeyFrame &KFc=*myMap->getKF(idKF);
					
					VectorXf dp(Dp.segment(7*k+1,6));
					if(verb_BA)std::cout<<"\t\tpose up "<<HomogeneousMatrix22(dp)<<std::endl;
					if(verb_BA)std::cout<<"\t\tscale up "<<scales_kfs[k]<<std::endl;
					KFc.setPose(pose_kfs[k]);
					for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
						KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
					
					HomogeneousMatrix22 relPose= KFc.getRelativePose();
					relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
					KFc.setRelativePose(relPose);
					
					HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
					relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
					KFc.setBestRelPose(relBestPose);
					
					scales_kfs[k]=1;
				}
			
				if(verb_BA)std::cout<<"\tupdate point positions"<<std::endl;
				//update position of points
				for(int k=0;k<FullWindow.size();k++)
				{
					int &idKF=FullWindow[k];
					//std::cout<<"\tidKF "<<idKF<<std::endl;
					KeyFrame &KFc=*myMap->getKF(idKF);
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
							KeyFrame &KFv=*myMap->getKF(kfView);
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
			if(verb_BA)std::cout<<"\tend loop iter"<<std::endl;
			
			/*{
				int idKFdisp=1;
				std::cout<<"Check BA points kf "<<idKFdisp<<std::endl;
				for(int p=0;p<myMap->getKF(idKFdisp)->getNbMapPoint();p++)
				{
					if(p<3)
					{
						MapPoint *point=myMap->getKF(idKFdisp)->getPtMapPoint(p);
						std::cout<<"point "<<p<<std::endl;
						std::cout<<"\tCoord : "<<point->getPosition().transpose()<<std::endl;
						std::cout<<"\tnbview "<<point->nbViews()<<std::endl;
						
						//check coord from local features
						for(int v=0;v<point->nbViews();v++)
						{
							int kfView=point->getView(v);
							int i1pView=point->getI1p(v);
							//get corresponding local feature:
							int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
							if(id_feat_v!=-1)
							{
								uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
								std::cout<<"\tCoord view : "<<feat_v.getLocalCoordinates().transpose()<<std::endl;
							}
						}
					}
				}
				//std::cout<<"Pose KFs "<<std::endl;
				//for(int k=0;k<nbOptimKf;k++)std::cout<<myMap->getKF(k)->getPose()<<std::endl;

			}*/
		}
		

		if(verb_BA)std::cout<<"\tLeave Optim"<<std::endl;
		//exit(1);
		
		/*{
			int idKFdisp=1;
			std::cout<<"Check BA points kf "<<idKFdisp<<std::endl;
			for(int p=0;p<myMap->getKF(idKFdisp)->getNbMapPoint();p++)
			{
				if(p<3)
				{
					MapPoint *point=myMap->getKF(idKFdisp)->getPtMapPoint(p);
					std::cout<<"point "<<p<<std::endl;
					std::cout<<"\tCoord : "<<point->getPosition().transpose()<<std::endl;
					std::cout<<"\tnbview "<<point->nbViews()<<std::endl;
					
					//check coord from local features
					for(int v=0;v<point->nbViews();v++)
					{
						int kfView=point->getView(v);
						int i1pView=point->getI1p(v);
						//get corresponding local feature:
						int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
						if(id_feat_v!=-1)
						{
							uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
							std::cout<<"\tCoord view : "<<feat_v.getLocalCoordinates().transpose()<<std::endl;
						}
					}
				}
			}
			//std::cout<<"Pose KFs "<<std::endl;
			//for(int k=0;k<nbOptimKf;k++)std::cout<<myMap->getKF(Inner(k))->getPose()<<std::endl;

		}*/
	}
}
	
}


struct dijkstra
{
	int kf_id;
	float rescale;
};

void MapOptimiser::optimiseScale(std::vector<int> &_innerWindowKFs,std::vector<scale_constraint> &mScaleConstraints,int nb_iter)
{
  	if(_innerWindowKFs.size()!=0)
	{
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	
	//scales of outerWindow needs to stay the same
	std::vector<float> outerWindowScales;for(int i=0;i<outerWindow.size();i++)outerWindowScales.push_back(1.);
	//set constraint as outerWindows with other scale
	for(int c=0;c<mScaleConstraints.size();c++)
	{
		for(int i=0;i<innerWindowKFs.size();i++)
		{
			if(mScaleConstraints[c].kf_id==innerWindowKFs[i])
			{
				outerWindow.push_back(innerWindowKFs[i]);
				outerWindowScales.push_back(mScaleConstraints[c].rescale);
				innerWindowKFs.erase(innerWindowKFs.begin()+i);
				break;
			}
		}
	}

	std::cout<<"innerWindowKFs = "<<std::endl;
	for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	std::cout<<std::endl;
	std::cout<<"outerWindow = "<<std::endl;
	for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	std::cout<<std::endl;
	
	if(outerWindow.size()==0 || innerWindowKFs.size()==0)
	{
		if(outerWindow.size()==0)
			coutRed<<"MapOptimiser::optimiseScale No fix frame around innerWindow found. Cannot proceed."<<endlRed;
	}
	else
	{
		int nbOptimKf=innerWindowKFs.size();
		bool verb_BA=true;
	
		//union of inner and outer
		std::vector<int> FullWindow=innerWindowKFs;
		for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
		
		//init LUT to get corresponding scale of KF in estimated vector scales_kfs
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		
		//get all the edges in FullWindow, will want scale to be constant on them
		std::vector<kf_edgeNew> list_edges=getEdgesInInnerWin(myMap,FullWindow);
		float min_weight_edge=list_edges[0].weight;
		for(int e=0;e<list_edges.size();e++)
		{
			kf_edgeNew &edge=list_edges[e];
			std::cout<<"edge : "<<edge.kf1<<" <=> "<<edge.kf2<<"  weight = "<<edge.weight<<std::endl;
			if(min_weight_edge>edge.weight)min_weight_edge=edge.weight;
		}
		
	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf+mScaleConstraints.size()];//+mScaleConstraints.size() will be used later to add constraints in scales_kfs
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		
		//TODO: first solution with scale = 1 except on constraint takes ages to propagate 
		//and stays stuck to local minimum => best would if we have two constraints only to
		//do max flow optimisation to find best cut=> one side = constraint one, other side = constraint two
		//then do optim.
		
		//or could compute shortest path to constraint and define init solution as constraint
		//distance could be defined as multiplication of inv edgeScore or as inv min edge score on path
		//mmm... quite tricky if (good path + scale) => crap path => (good path no scale) => crap path => (good path and scale)
		//what s the most likely scale of middle good path ?
		//if distance defined by inv min edge score on path then will init with side with less crapy path...
		//would be good
		
		//first scale
		float ScalesInit[nbOptimKf];
		for(int i=0;i<nbOptimKf;i++)ScalesInit[i]=-1;//+1 to be sure its bigger
		
		{
			float dist_min[nbOptimKf];
			//init dist min with dist max = inv min edge
			float max_dist=1./min_weight_edge+1;
			for(int i=0;i<nbOptimKf;i++)dist_min[i]=max_dist;//+1 to be sure its bigger
			
			
			//start from each constraint outerWindow with scales outerWindowScales and propagate scale on neigbors if distance is inferior to current dist_min
			for(int constr=0;constr<outerWindow.size();constr++)
			{
				//define source kf from which to explore neigbors until everything is explored
				int kf_start=outerWindow[constr];
				float scale_start=outerWindowScales[constr];
				
				std::vector<int> currentToCheckNeigbor;
				std::vector<int> toBeAddedInCurrentToCheckNeigbor;
				currentToCheckNeigbor.push_back(kf_start);
				//std::vector<int> toCheckKf=innerWindowKFs;

				std::vector<int> alreadyChecked;
				std::vector<int> toBeAddedInalreadyChecked;

				do
				{
					for(int i=0;i<currentToCheckNeigbor.size();i++)
					{
						int &idi=currentToCheckNeigbor[i];
						KeyFrame &kfi=*myMap->getKF(idi);
						
						int idOpt=LUTidKFtoInnerWindow[idi];
						float dist_current;
						if(idOpt==-1)//outside window opt=> is source
							dist_current=0;
						else
							dist_current=dist_min[idOpt];
						
						//go through all neigbors
						for(int n=0;n<kfi.getNbNeigbours();n++)
						{
							NeigbourKFNew &neigbor=*kfi.getPtNeigbour(n);
							//if hasn t been already checked
							if(std::find(alreadyChecked.begin(), alreadyChecked.end(), neigbor.neighboring_kf)==alreadyChecked.end())
							{
								//define new dist as min distance current and inv edge score
								float inv_score_edge=1./neigbor.edgeScore;
								float dist_neigbor_from_current=(inv_score_edge>dist_current)?inv_score_edge:dist_current;
								
								//if new distance of neigbor smaller than distance that was stored in it then we found 
								//that this constraint is closer=> change
								int idOptNeigb=LUTidKFtoInnerWindow[neigbor.neighboring_kf];
								if(idOptNeigb!=-1)
								  if(dist_min[idOptNeigb]>dist_neigbor_from_current)
								  {
								    dist_min[idOptNeigb]=dist_neigbor_from_current;
								    ScalesInit[idOptNeigb]=scale_start;
								    toBeAddedInCurrentToCheckNeigbor.push_back(neigbor.neighboring_kf);
								  }
								
							}
							toBeAddedInalreadyChecked.push_back(neigbor.neighboring_kf);
							
						}
					}
					
					//update kf to be checked on last step
					currentToCheckNeigbor.clear();
					for(int i=0;i<toBeAddedInCurrentToCheckNeigbor.size();i++)currentToCheckNeigbor.push_back(toBeAddedInCurrentToCheckNeigbor[i]);
					toBeAddedInCurrentToCheckNeigbor.clear();
					//update kf to be checked on last step
					for(int i=0;i<toBeAddedInalreadyChecked.size();i++)
						if(std::find(alreadyChecked.begin(), alreadyChecked.end(), toBeAddedInalreadyChecked[i])==alreadyChecked.end())
						  alreadyChecked.push_back(toBeAddedInalreadyChecked[i]);
					toBeAddedInalreadyChecked.clear();
					  
					
					
				}
				while(currentToCheckNeigbor.size()!=0);
			}
		}
		//init scales with closest fix scale
		for(int k=0;k<nbOptimKf;k++)
			scales_kfs[k]=ScalesInit[k];
			
		
		  
		/*std::cout<<"scales = "<<std::endl;
		for(int i=0;i<nbOptimKf;i++)
			std::cout<<scales_kfs[i] <<" ";
		std::cout<<std::endl<<std::endl;*/
		
		
		//get through all neigbors of outerWindowInit, everytime 
				
		int LUTidKFtoOuterWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoOuterWindow[i]=-1;
		for(int i=0;i<outerWindow.size();i++)LUTidKFtoOuterWindow[outerWindow[i]]=i;
		
		for(int iter=0;iter<nb_iter;iter++)
		{
		  
			int nb_params=nbOptimKf;
			VectorXf Jte(nb_params);Jte.setZero();
			//MatrixXf H(nb_params,nb_params);H.setZero();
			
			float residue=0;
			
			//create error from edges
			for(int e=0;e<list_edges.size();e++)
			{
				kf_edgeNew &edge=list_edges[e];

				int id_opt1=LUTidKFtoInnerWindow[edge.kf1];
				int id_opt2=LUTidKFtoInnerWindow[edge.kf2];
				
				//get scale current of two kf from that
				float scale_c1;
				if(id_opt1!=-1)scale_c1=scales_kfs[id_opt1];//being optimised=>get from estimated scales
				else scale_c1=outerWindowScales[LUTidKFtoOuterWindow[edge.kf1]];//else set as fixed in outerWindowScales
				
				float scale_c2;
				if(id_opt2!=-1)scale_c2=scales_kfs[id_opt2];//being optimised=>get from estimated scales
				else scale_c2=outerWindowScales[LUTidKFtoOuterWindow[edge.kf2]];
				
				//scale diff, want it to be equal to 0, with confidence = weight edge
				float scale_diff=scale_c1-scale_c2;
				residue+=scale_diff*scale_diff;
				
				//update jacobians
				if(id_opt1!=-1)
					Jte[id_opt1]+=edge.weight*scale_diff;
				if(id_opt2!=-1)
					Jte[id_opt2]-=edge.weight*scale_diff;
				
				//update Hessian
				/*if(id_opt1!=-1)
					H(id_opt1,id_opt1)+=edge.weight;
				if(id_opt2!=-1)
					H(id_opt2,id_opt2)+=edge.weight;
				if(id_opt1!=-1 && id_opt2!=-1)
				{
					H(id_opt1,id_opt2)+=edge.weight;
					H(id_opt2,id_opt1)+=edge.weight;				  
				}*/
			}
			
			/*if(iter==0)
			{
			std::cout<<"H = "<<std::endl;
			std::cout<<H<<std::endl;
			}*/
			
			//std::cout<<"Jac = "<<std::endl;
			//std::cout<<Jte.transpose()<<std::endl;
			
			//std::cout<<"residue = "<<residue<<std::endl;
			
			//do update
			/*Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
			  
				VectorXf Ds(nb_params);
				Ds=-0.3*(invH*Jte);
				
				for(int i=0;i<nb_params;i++)
					scales_kfs[i]+=Ds[i];
			}*/
			
			//steppest gradient method
			VectorXf Ds(nb_params);
			Ds=-1e-6*Jte;
			
			for(int i=0;i<nb_params;i++)
				scales_kfs[i]+=Ds[i];			
			
			
	
		}
		
		//add constraints:
		for(int c=0;c<mScaleConstraints.size();c++)
		{
			scales_kfs[nbOptimKf+c]=mScaleConstraints[c].rescale;
			innerWindowKFs.push_back(mScaleConstraints[c].kf_id);
		}
		//update LUT
		nbOptimKf=nbOptimKf+mScaleConstraints.size();
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;

		
		std::cout<<"scales = "<<std::endl;
		for(int i=0;i<nbOptimKf;i++)
			std::cout<<scales_kfs[i] <<" ";
		std::cout<<std::endl<<std::endl;
			  
		
		
		//end iter, have to change local features in KF
		std::cout<<"rescale local features"<<std::endl;
		for(int k=0;k<nbOptimKf;k++)
		{
			int &idKF=innerWindowKFs[k];
			if(verb_BA)std::cout<<"\t\tupdate "<<idKF<<std::endl;
			KeyFrame &KFc=*myMap->getKF(idKF);
			
			for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
				KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
		}
				
		std::cout<<"rescale rel poses"<<std::endl;
		//have to rescale kf relative positions
		for(int e=0;e<list_edges.size();e++)
		{
			kf_edgeNew &edge=list_edges[e];

			int id_opt1=LUTidKFtoInnerWindow[edge.kf1];
			int id_opt2=LUTidKFtoInnerWindow[edge.kf2];
			
			//get scale current of two kf from that
			float scale_c1;
			if(id_opt1!=-1)scale_c1=scales_kfs[id_opt1];//being optimised=>get from estimated scales
			else scale_c1=outerWindowScales[LUTidKFtoOuterWindow[edge.kf1]];//else set as fixed in outerWindowScales
			
			float scale_c2;
			if(id_opt2!=-1)scale_c2=scales_kfs[id_opt2];//being optimised=>get from estimated scales
			else scale_c2=outerWindowScales[LUTidKFtoOuterWindow[edge.kf2]];
			
			//rescale edge by mean scales
			float meanScale=0.5*(scale_c1+scale_c2);
			
			KeyFrame &kf1=*myMap->getKF(edge.kf1);
			for(int n=0;n<kf1.getNbNeigbours();n++)
			{
				NeigbourKFNew &neigbor=*kf1.getPtNeigbour(n);
				if(neigbor.neighboring_kf==edge.kf2)
					neigbor.relative_poses.set_translation(meanScale*neigbor.relative_poses.get_translation());
			}
			
			KeyFrame &kf2=*myMap->getKF(edge.kf2);
			for(int n=0;n<kf2.getNbNeigbours();n++)
			{
				NeigbourKFNew &neigbor=*kf2.getPtNeigbour(n);
				if(neigbor.neighboring_kf==edge.kf1)
					neigbor.relative_poses.set_translation(meanScale*neigbor.relative_poses.get_translation());
			}
			
		}
		
		std::cout<<"rescale current rel poses and best stereo pair pose"<<std::endl;
		//have to rescale kf relative current positions
		for(int k=0;k<nbOptimKf;k++)
		{
			int &idKF=innerWindowKFs[k];
			KeyFrame &KFc=*myMap->getKF(idKF);
		 	HomogeneousMatrix22 relPose= KFc.getRelativePose();
			relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
			KFc.setRelativePose(relPose);
			
		 	HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
			relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
			KFc.setBestRelPose(relBestPose);
		}
	
		
		std::cout<<"update kf positions"<<std::endl;
		//have to update kf positions
		//tricky one again, here we decide to estimate position of KF using path
		//(ie composition) of relative poses and one fixed KF pose
		//reset inner and outer windows to original thing (constraint on scale doesn t make kf fixed)
		outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
		
		//start from all outerWindow do iterative process to explore neigbors and check if current checked kf is closer to current fixed 
		{
			float dist_min[nbOptimKf];
			//init dist min with dist max = inv min edge
			float max_dist=1./min_weight_edge+1;
			for(int i=0;i<nbOptimKf;i++)dist_min[i]=max_dist;//+1 to be sure its bigger
			
			
			//start from each constraint outerWindow with scales outerWindowScales and propagate scale on neigbors if distance is inferior to current dist_min
			for(int o=0;o<outerWindow.size();o++)
			{
				//define source kf from which to explore neigbors until everything is explored
				int kf_start=outerWindow[o];
				
				std::vector<int> currentToCheckNeigbor;
				std::vector<int> toBeAddedInCurrentToCheckNeigbor;
				currentToCheckNeigbor.push_back(kf_start);
				//std::vector<int> toCheckKf=innerWindowKFs;

				std::vector<int> alreadyChecked;
				std::vector<int> toBeAddedInalreadyChecked;

				do
				{
					for(int i=0;i<currentToCheckNeigbor.size();i++)
					{
						int &idi=currentToCheckNeigbor[i];
						KeyFrame &kfi=*myMap->getKF(idi);
						
						int idOpt=LUTidKFtoInnerWindow[idi];
						float dist_current;
						if(idOpt==-1)//outside window opt=> is source
							dist_current=0;
						else
							dist_current=dist_min[idOpt];
						
						//go through all neigbors
						for(int n=0;n<kfi.getNbNeigbours();n++)
						{
							NeigbourKFNew &neigbor=*kfi.getPtNeigbour(n);
							//if hasn t been already checked
							if(std::find(alreadyChecked.begin(), alreadyChecked.end(), neigbor.neighboring_kf)==alreadyChecked.end())
							{
								//define new dist as min distance current and inv edge score
								float inv_score_edge=1./neigbor.edgeScore;
								float dist_neigbor_from_current=inv_score_edge+dist_current;
								
								//if new distance of neigbor smaller than distance that was stored in it then we found 
								//that this constraint is closer=> change
								int idOptNeigb=LUTidKFtoInnerWindow[neigbor.neighboring_kf];
								if(idOptNeigb!=-1)
								  if(dist_min[idOptNeigb]>dist_neigbor_from_current)
								  {
								    dist_min[idOptNeigb]=dist_neigbor_from_current;
								   // ScalesInit[idOptNeigb]=scale_start;
								    myMap->getKF(neigbor.neighboring_kf)->setPose(neigbor.relative_poses.inverse()*kfi.getPose());
								    toBeAddedInCurrentToCheckNeigbor.push_back(neigbor.neighboring_kf);
								  }
								
							}
							toBeAddedInalreadyChecked.push_back(neigbor.neighboring_kf);
							
						}
					}
					
					//update kf to be checked on last step
					currentToCheckNeigbor.clear();
					for(int i=0;i<toBeAddedInCurrentToCheckNeigbor.size();i++)currentToCheckNeigbor.push_back(toBeAddedInCurrentToCheckNeigbor[i]);
					toBeAddedInCurrentToCheckNeigbor.clear();
					//update kf to be checked on last step
					for(int i=0;i<toBeAddedInalreadyChecked.size();i++)
						if(std::find(alreadyChecked.begin(), alreadyChecked.end(), toBeAddedInalreadyChecked[i])==alreadyChecked.end())
						  alreadyChecked.push_back(toBeAddedInalreadyChecked[i]);
					toBeAddedInalreadyChecked.clear();
					  
					
					
				}
				while(currentToCheckNeigbor.size()!=0);
			}
		}
		
		
		//have to reestimate map point posisions
		optimiseMapPoints(FullWindow);
		
	}
	}
}

void MapOptimiser::optimiseMapPoints(std::vector<int> &_innerWindowKFs)
{
	std::cout<<"optimiseMapPoints"<<std::endl;
	for(int k=0;k<_innerWindowKFs.size();k++)
	{
		int &idKF=_innerWindowKFs[k];
		KeyFrame &KFc=*myMap->getKF(idKF);
		for(int p=0;p<KFc.getNbMapPoint();p++)
		{
			MapPoint *point=KFc.getPtMapPoint(p);
			if(point->isUsed())//if point is good
			{
				Vector3f newPos;newPos.setZero();
				int nbViews=0;
				for(int v=0;v<point->nbViews();v++)
				{
					int kfView=point->getView(v);
					int i1pView=point->getI1p(v);
					//get corresponding local feature:
					int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
					if(id_feat_v!=-1)
					{
						uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
						//compute error position point versus position local feature scaled 
						Vector3f pointFromLocalFeat=myMap->getKF(kfView)->getPose().inverse()*feat_v.getLocalCoordinates();
						newPos+=pointFromLocalFeat;
						nbViews++;
					}
				}
				point->updatePosition(newPos/nbViews);
			}
		}				
	}
  
}


/*
void MapOptimiser::optimiseScale(std::vector<int> &_innerWindowKFs,std::vector<scale_constraint> &mScaleConstraints,int nb_iter)
{
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);

	std::cout<<"innerWindowKFs = "<<std::endl;
	for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;
	std::cout<<"outerWindow = "<<std::endl;
	for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  "<<std::endl;
	
	if(outerWindow.size()==0 || innerWindowKFs.size()==0)
	{
		if(outerWindow.size()==0)
			coutRed<<"MapOptimiser::optimiseScale No fix frame around innerWindow found. Cannot proceed."<<endlRed;
	}
	else
	{
		int nbOptimKf=innerWindowKFs.size();
		bool verb_BA=true;
	
		//union of inner and outer
		std::vector<int> FullWindow=innerWindowKFs;
		for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
		
		//get all the edges in FullWindow, will want scale to be constant on them
		std::vector<kf_edgeNew> list_edges=getEdgesInInnerWin(FullWindow);
		for(int e=0;e<list_edges.size();e++)
		{
			kf_edgeNew &edge=list_edges[e];
			std::cout<<"edge : "<<edge.kf1<<" <=> "<<edge.kf2<<"  weight = "<<edge.weight<<std::endl;
		}
		
	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		
		//init constrained scales with their desired values
		for(int c=0;c<mScaleConstraints.size();c++)scales_kfs[mScaleConstraints[c].kf_id]=mScaleConstraints[c].rescale;
		
		//init LUT to get corresponding scale of KF in estimated vector scales_kfs
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		
		for(int iter=0;iter<nb_iter;iter++)
		{
		  
			int nb_params=nbOptimKf;
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			float residue=0;
			
			//create error from edges
			for(int e=0;e<list_edges.size();e++)
			{
				kf_edgeNew &edge=list_edges[e];

				int id_opt1=LUTidKFtoInnerWindow[edge.kf1];
				int id_opt2=LUTidKFtoInnerWindow[edge.kf2];
				
				//get scale current of two kf from that
				float scale_c1=1.;//if KF is fixed then scale stays = 1
				if(id_opt1!=-1)scale_c1=scales_kfs[id_opt1];//else it is what is being optimised
				float scale_c2=1.;//if KF is fixed then scale stays = 1
				if(id_opt2!=-1)scale_c2=scales_kfs[id_opt2];//else it is what is being optimised
				
				//scale diff, want it to be equal to 0, with confidence = weight edge
				float scale_diff=scale_c1-scale_c2;
				residue+=scale_diff*scale_diff;
				
				//update jacobians
				if(id_opt1!=-1)
					Jte[id_opt1]+=edge.weight*scale_diff;
				if(id_opt2!=-1)
					Jte[id_opt2]-=edge.weight*scale_diff;
				
				//update Hessian
				if(id_opt1!=-1)
					H(id_opt1,id_opt1)+=edge.weight;
				if(id_opt2!=-1)
					H(id_opt2,id_opt2)+=edge.weight;
				if(id_opt1!=-1 && id_opt2!=-1)
				{
					H(id_opt1,id_opt2)+=edge.weight;
					H(id_opt2,id_opt1)+=edge.weight;				  
				}
			}
			
			//create error from constraints
			float Lambda=1e5;//make it very important
			for(int c=0;c<mScaleConstraints.size();c++)
			{
				scale_constraint &constraint=mScaleConstraints[c];
				
				int id_opt=LUTidKFtoInnerWindow[constraint.kf_id];
				if(id_opt!=-1)
				{
					float scale_c=scales_kfs[id_opt];//else it is what is being optimised
					float scale_diff=scale_c-constraint.rescale;
					
					//update jacobian and Hessian
					Jte[id_opt]+=Lambda*scale_diff;
				  
					//update jacobian and Hessian
					H(id_opt,id_opt)+=Lambda;
					
					residue+=scale_diff*scale_diff;
				}
			}
			
			if(iter==0)
			{
			std::cout<<"H = "<<std::endl;
			std::cout<<H<<std::endl;
			}
			
			std::cout<<"Jac = "<<std::endl;
			std::cout<<Jte.transpose()<<std::endl;
			
			std::cout<<"residue = "<<residue<<std::endl;
			
			//do update
			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
			  
				VectorXf Ds(nb_params);
				Ds=-0.3*(invH*Jte);
				
				for(int i=0;i<nb_params;i++)
					scales_kfs[i]+=Ds[i];
			}
			
			//steppest gradient method
			//VectorXf Ds(nb_params);
			//Ds=-1e-6*Jte;
			
			//for(int i=0;i<nb_params;i++)
			//	scales_kfs[i]+=Ds[i];			
			
		
			std::cout<<"scales = "<<std::endl;
			for(int i=0;i<nb_params;i++)
				std::cout<<scales_kfs[i] <<" ";
			std::cout<<std::endl<<std::endl;
			  
	
		}
		
		
		//end iter, have to change local features in KF
	
		//have to rescale kf relative positions
		
		//have to update kf positions
		
		//have to reestimate map point posisions
		
	}
}*/
void MapOptimiser::optimiseInnerWindowRobust2(std::vector<int> &_innerWindowKFs,int nb_iter)
{
	if(_innerWindowKFs.size()!=0)
	{
	coutGreen<<"###############################################"<<endlGreen;
	coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	//std::cout<<"outerWindow = "<<std::endl;
	//for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  "<<std::endl;
	//need one KF to be fixed
	//if empty then need to fix one frame, take first one out
	if(outerWindow.size()==0)
	{
		outerWindow.push_back(*innerWindowKFs.begin());
		innerWindowKFs.erase(innerWindowKFs.begin());
	}

	//union of inner and outer
	std::vector<int> FullWindow=innerWindowKFs;
	for(int k=0;k<outerWindow.size();k++)FullWindow.push_back(outerWindow[k]);
	
	//check windows:
	//std::cout<<"FullWindow = "<<std::endl;
	//for(int i=0;i<FullWindow.size();i++)std::cout<<FullWindow[i]<<"  "<<std::endl;
	std::cout<<"innerWindowKFs = "<<std::endl;
	for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	std::cout<<std::endl;
	std::cout<<"outerWindow = "<<std::endl;
	for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	std::cout<<std::endl<<std::endl;
	
	bool verb_BA=false;
	int nbOptimKf=innerWindowKFs.size();
	
		
	if(nbOptimKf>0)
	{

	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		float scales_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
		HomogeneousMatrix22 pose_kfs[nbOptimKf];
		for(int k=0;k<nbOptimKf;k++)pose_kfs[k]=myMap->getKF(innerWindowKFs[k])->getPose();
		
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		  
		float mLMLambda = 0.0001;//before LevMarConstant
		float mdLambdaFactor = 2.0;
		
		if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
		for(int iter=0;iter<nb_iter;iter++)
		{
			if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			//need to put in minimisation all the points linked by features of innerWindow and outerWindow
			//that should be all the points defined in this keyframes
			
			//get TukeyScalar
			float sigma_tukey[FullWindow.size()];
			for(int k=0;k<FullWindow.size();k++)
			{
				std::vector<float> vdErrorSquared;
				int &idKF=FullWindow[k];
				KeyFrame &KFc=*myMap->getKF(idKF);
				for(int p=0;p<KFc.getNbMapPoint();p++)
				{
					//now we have the point and each of its views
					//each view has to be used to compute the sum to minimise
					MapPoint *point=KFc.getPtMapPoint(p);
					if(point->isUsed())//if point is good
					{
						for(int v=0;v<point->nbViews();v++)
						{
							int kfView=point->getView(v);
							int i1pView=point->getI1p(v);
							//get corresponding local feature:
							int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
							if(id_feat_v!=-1)
							{
								uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
								int id_opt=LUTidKFtoInnerWindow[kfView];
								if(id_opt!=-1)//if corresponding keyframe is to be optimised
								{
									//compute error position point versus position local feature scaled 
									Vector3f pointInCurrentKF=pose_kfs[id_opt]*point->getPosition();
									Vector3f pointFromLocalFeat=scales_kfs[id_opt]*feat_v.getLocalCoordinates();
									Vector3f positionError= pose_kfs[id_opt]*point->getPosition()-pointFromLocalFeat;
									vdErrorSquared.push_back(positionError.transpose()*positionError);
								}
							}
						}
					}
				}
				if(vdErrorSquared.size()>0)
					sigma_tukey[k]=getSigmaSquared(vdErrorSquared);
				else	
					sigma_tukey[k]=0;
			}
		
			//residue just to check if error goes down
			std::vector<kf_loc_jacobian> list_jacobian;
			float residue=0;
			float totalWeight=0;
			for(int k=0;k<FullWindow.size();k++)
			{
				int &idKF=FullWindow[k];
				KeyFrame &KFc=*myMap->getKF(idKF);
				if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb points : "<<KFc.getNbMapPoint()<<std::endl;
				for(int p=0;p<KFc.getNbMapPoint();p++)
				{
					
					//now we have the point and each of its views
					//each view has to be used to compute the sum to minimise
					MapPoint *point=KFc.getPtMapPoint(p);
					if(point->isUsed())//if point is good
					{
						for(int v=0;v<point->nbViews();v++)
						{
							int kfView=point->getView(v);
							int i1pView=point->getI1p(v);
							//get corresponding local feature:
							int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
							if(id_feat_v!=-1)
							{
								uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
								int id_opt=LUTidKFtoInnerWindow[kfView];
								if(id_opt!=-1)//if corresponding keyframe is to be optimised
								{
									//compute error position point versus position local feature scaled 
									Vector3f pointInCurrentKF=pose_kfs[id_opt]*point->getPosition();
									Vector3f pointFromLocalFeat=scales_kfs[id_opt]*feat_v.getLocalCoordinates();
									Vector3f positionError= pose_kfs[id_opt]*point->getPosition()-pointFromLocalFeat;
									
									float TukeyCoef=squareRootTukey(positionError.squaredNorm(),sigma_tukey[k]);
									if(TukeyCoef==0 && iter>2)//give two iterations to have acceptable error, after that consider as outlier and will be suppressed
										feat_v.nb_outlier++;
									
									if(TukeyCoef>0 && feat_v.scoreFundamentalOrigin>0)
									{
										
										residue+=positionError.squaredNorm();
										
										//create jacobian
										kf_loc_jacobian fJacobian;
										fJacobian.pos_error=positionError;
										fJacobian.opt_id=id_opt;
										
										//compute jacobians:
										//jacobian wrt kf position: consider that search for update dp to compose with pose_kfs[id_opt] to minise equation
										
										//has to be 3by6: diff of pos error 3 by 6 dof
										MatrixXf de_dp_t(3,6);
										de_dp_t.block(0,0,3,3)=Matrix3f::Identity();
										de_dp_t.block(0,3,3,3)=GetSkew(-pointInCurrentKF);
																		
										fJacobian.de_dp=de_dp_t;
										
										
										//jacobian wrt kf scale
										//Xlocal2=(1+ds)s*Xlocal=>
										fJacobian.de_ds=-scales_kfs[id_opt]*feat_v.getLocalCoordinates();
										
										
										
										//can use weight for this measure as reconstruction angle point or fundamental matrix score of origin
										//fJacobian.weight=TukeyCoef*feat_v.scoreFundamentalOrigin;
										fJacobian.weight=TukeyCoef*feat_v.recAngle;
										totalWeight+=TukeyCoef*feat_v.recAngle;
										
										list_jacobian.push_back(fJacobian);
									}
								}
							}
							
						}
					}
				}				
			}
			if(verb_BA)coutRed<<"\tresitue BA = "<<residue<<endlRed;
			if(verb_BA)std::cout<<"\tresitue BA normalized= "<<residue/list_jacobian.size()<<std::endl;
			
			//compute update from jacobians
			int nb_params=7*nbOptimKf;
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			if(verb_BA)std::cout<<"\taccumulate in J and H"<<std::endl;
			for(int j=0;j<list_jacobian.size();j++)
			{
				kf_loc_jacobian &fJacobian=list_jacobian[j];
				
				//scale first
				Jte[7*fJacobian.opt_id]+=fJacobian.weight * fJacobian.de_ds.transpose()*fJacobian.pos_error;
				//then kf position
				Jte.segment(7*fJacobian.opt_id+1,6)+=fJacobian.weight * fJacobian.de_dp.transpose()*fJacobian.pos_error;
				
				//Hessian:
				H(7*fJacobian.opt_id,7*fJacobian.opt_id)+=fJacobian.weight * fJacobian.de_ds.transpose() * fJacobian.de_ds;				
				H.block(7*fJacobian.opt_id+1,7*fJacobian.opt_id+1,6,6)+=fJacobian.weight * fJacobian.de_dp.transpose() * fJacobian.de_dp;				
			}
			
			if(verb_BA)std::cout<<"\tcompute update"<<std::endl;
			//std::cout<<"\tJte = "<<std::endl;
			//std::cout<<Jte<<std::endl;
			//std::cout<<"\tH = "<<std::endl;
			//std::cout<<H<<std::endl;
			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			if(isnan(invH(0,0)))
			{
				std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
			  
				VectorXf Dp(nb_params);
				Dp=-gain*(invH*Jte);
				
				if(verb_BA)std::cout<<"\tapply update"<<std::endl;
				//apply update to kf
				for(int j=0;j<nbOptimKf;j++)
				{
					if(verb_BA)std::cout<<"\tapply update scale"<<std::endl;
					float ds=Dp[7*j];
					scales_kfs[j]=(1.+ds)*scales_kfs[j];
					
					if(verb_BA)std::cout<<"\tapply update pose"<<std::endl;
					VectorXf dp(Dp.segment(7*j+1,6));
					pose_kfs[j]=HomogeneousMatrix22(dp)* pose_kfs[j];
				}
				
				float residue_after=0;
				float totalWeight_after=0;
				for(int k=0;k<FullWindow.size();k++)
				{
					int &idKF=FullWindow[k];
					KeyFrame &KFc=*myMap->getKF(idKF);
					if(verb_BA)std::cout<<"\tidKF "<<idKF<<" Nb points : "<<KFc.getNbMapPoint()<<std::endl;
					for(int p=0;p<KFc.getNbMapPoint();p++)
					{
						
						//now we have the point and each of its views
						//each view has to be used to compute the sum to minimise
						MapPoint *point=KFc.getPtMapPoint(p);
						if(point->isUsed())//if point is good
						{
							for(int v=0;v<point->nbViews();v++)
							{
								int kfView=point->getView(v);
								int i1pView=point->getI1p(v);
								//get corresponding local feature:
								int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
								if(id_feat_v!=-1)
								{
									uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
									int id_opt=LUTidKFtoInnerWindow[kfView];
									if(id_opt!=-1)//if corresponding keyframe is to be optimised
									{
										//compute error position point versus position local feature scaled 
										Vector3f pointInCurrentKF=pose_kfs[id_opt]*point->getPosition();
										Vector3f pointFromLocalFeat=scales_kfs[id_opt]*feat_v.getLocalCoordinates();
										Vector3f positionError= pose_kfs[id_opt]*point->getPosition()-pointFromLocalFeat;
										
										float TukeyCoef=squareRootTukey(positionError.squaredNorm(),sigma_tukey[k]);
										if(TukeyCoef==0 && iter>2)//give two iterations to have acceptable error, after that consider as outlier and will be suppressed
											feat_v.nb_outlier++;
										
										if(TukeyCoef>0 && feat_v.scoreFundamentalOrigin>0)
										{
											residue_after+=positionError.squaredNorm();
											totalWeight_after+=TukeyCoef*feat_v.recAngle;
											
										}
									}
								}
								
							}
						}
					}				
				}
				
				if(residue_after/totalWeight_after<residue/totalWeight)
				{
					mdLambdaFactor = 2.0;
					mLMLambda *= 0.3;
				
					if(verb_BA)std::cout<<"\tupdate map with new estimates"<<std::endl;
					//update KFs with new scale and pose //make it faster to compute point position of do that here in loop
					for(int k=0;k<nbOptimKf;k++)
					{
						if(verb_BA)std::cout<<"\t\tupdate "<<k<<std::endl;
						int &idKF=innerWindowKFs[k];
						KeyFrame &KFc=*myMap->getKF(idKF);
						
						VectorXf dp(Dp.segment(7*k+1,6));
						if(verb_BA)std::cout<<"\t\tpose up "<<HomogeneousMatrix22(dp)<<std::endl;
						if(verb_BA)std::cout<<"\t\tscale up "<<scales_kfs[k]<<std::endl;
						KFc.setPose(pose_kfs[k]);
						for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
							KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
						
						
						HomogeneousMatrix22 relPose= KFc.getRelativePose();
						relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
						KFc.setRelativePose(relPose);
						
						HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
						relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
						KFc.setBestRelPose(relBestPose);
						
						scales_kfs[k]=1;
					}
				
					if(verb_BA)std::cout<<"\tupdate point positions"<<std::endl;
					//update position of points
					for(int k=0;k<FullWindow.size();k++)
					{
						int &idKF=FullWindow[k];
						//std::cout<<"\tidKF "<<idKF<<std::endl;
						KeyFrame &KFc=*myMap->getKF(idKF);
						for(int p=0;p<KFc.getNbMapPoint();p++)
						{
							
							//if(p<10)std::cout<<"\t\tpoint "<<p<<std::endl;
							//now we have the point and each of its views
							//each view has to be used to compute the sum to minimise
							MapPoint *point=KFc.getPtMapPoint(p);
							
							float newWeight=0;
							Vector3f newPosition=Vector3f(0,0,0);
							if(point->isUsed())//if point is good
							{
								for(int v=0;v<point->nbViews();v++)
								{
									//if(p<10)std::cout<<"\t\t\tview = "<<v<<std::endl;
									int kfView=point->getView(v);
									//if(p<10)std::cout<<"\t\t\tkfview = "<<kfView<<std::endl;
									int i1pView=point->getI1p(v);
									KeyFrame &KFv=*myMap->getKF(kfView);
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
				}
				else
				{
				  	mLMLambda = mLMLambda * mdLambdaFactor;
					mdLambdaFactor = mdLambdaFactor * 2;
					for(int k=0;k<nbOptimKf;k++)
					{
						int &idKF=innerWindowKFs[k];
						KeyFrame &KFc=*myMap->getKF(idKF);
						pose_kfs[k]=KFc.getPose();
						scales_kfs[k]=1;
					}
				}
			}
			if(verb_BA)std::cout<<"\tend loop iter"<<std::endl;
			
		}
		
		//remove outliers
		for(int k=0;k<FullWindow.size();k++)
		{
			int &idKF=FullWindow[k];
			KeyFrame &KFc=*myMap->getKF(idKF);
			for(int p=0;p<KFc.getNbMapPoint();p++)
			{
				MapPoint *point=KFc.getPtMapPoint(p);
				if(point->isUsed())//if point is good
				{
					for(int v=0;v<point->nbViews();v++)
					{
						int kfView=point->getView(v);
						int i1pView=point->getI1p(v);
						//get corresponding local feature:
						int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
						if(id_feat_v!=-1)
						{
							uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
							if(feat_v.nb_outlier!=0)
							{
								//remove view
								point->removeViewId(v);
								//remove match from feature
								feat_v.matched=false;
								v--;
							}
						}
					}
					//check if point has more than one view, if not remove it
					if(point->nbViews()<=1)
					{
						if(point->nbViews()==1)
						{
							int kfView=point->getView(0);
							int i1pView=point->getI1p(0);
							int id_feat_v=myMap->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
							uptoscaleFeature &feat_v=*myMap->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
							feat_v.matched=false;
						}
						point->removeViews();
						point->setAsBad();
					}
				}
			}
		}
		
		

		if(verb_BA)std::cout<<"\tLeave Optim"<<std::endl;
		//exit(1);

	}
	}
	
}
void MapOptimiser::getRelativePoseAndScale(int _kf1,int _kf2,float &optRelScale,HomogeneousMatrix22 &optRelPose,float &infoScale,MatrixXf &infoPose)
{
	bool verb_BA=false;
	
	Camera *myCam=myMap->getCamera();
		
	KeyFrame &KFc=*myMap->getKF(_kf1);
	KeyFrame &KFn=*myMap->getKF(_kf2);

	if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
	//we are going to optimise over scale of each KF and relative pose between KF
	//only the pose and scale of second kf are optimised
	float scales_kf2=1.;
	HomogeneousMatrix22 pose_kf2=KFn.getPose();
	int nbOptimKf=1;
			
	if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
	int nb_iter=50;
	for(int iter=0;iter<nb_iter;iter++)
	{
		//if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
		//if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
		std::cout<<"\t###############################################"<<std::endl;
		std::cout<<"\tIteration "<<iter<<std::endl;
		//need to put in minimisation all the points linked by features of innerWindow and outerWindow
		//that should be all the points defined in this keyframes
		std::vector<kf_loc_jacobian2> list_jacobian_scale;
		
		HomogeneousMatrix22 pose_c=KFc.getPose();
		HomogeneousMatrix22 pose_n=pose_kf2;		
		HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
		
		//to get from relPose to deriv wrt pose_c and pose_n
		MatrixXf M1(6,6);
		VectorXf logRelPose=relPose.get_p();
		Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logRelPose[j];
		Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logRelPose[j+3];
		M1.block(0,0,3,3)=-GetSkew(Dw);
		M1.block(3,3,3,3)=-GetSkew(Dw);
		M1.block(0,3,3,3)=-GetSkew(Dt);
		M1.block(3,0,3,3).setZero();
		MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1);
		
		float residue=0;
		float nb_residue=0;

		//get neigbors
		for(int n=0;n<KFc.getNbNeigbours();n++)
		{
			NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
			if(neigbor.neighboring_kf==_kf2)//can do optim only if kf2 is neighbor of kf1
			{
				
				if(verb_BA)std::cout<<"\tnb matches = "<<neigbor.matches.size()<<std::endl;
				for(int m=0;m<neigbor.matches.size();m++)
				{
					p_match &mMatch=neigbor.matches[m];
					//if matches are linked to local features then can use them for pose and scale alignment
					int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
					int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
					if(idFeat_c!=-1 && idFeat_n!=-1)
					{
						//get current estimated scales
							
						float scale_n=scales_kf2;
							
						uptoscaleFeature *feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
						uptoscaleFeature *feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
						Vector3f scaleError=feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates());
						residue+=scaleError.squaredNorm();
						nb_residue++;
						
						
						//deriv with respect to relPose:
						MatrixXf jac_dp(3,6);
						jac_dp.block(0,0,3,3)=-Matrix3f::Identity();
						jac_dp.block(0,3,3,3)=GetSkew(relPose*(scale_n*feat_n->getLocalCoordinates()));
					
						
						kf_loc_jacobian2 mJacScale;
						mJacScale.pos_error=scaleError;
						mJacScale.weight=1.;
						
						
						mJacScale.opt_idn=0;
						//mJacScale.de_dpn=jac_dp*dHErr_dpn;
						mJacScale.de_dpn=jac_dp;//instead will multiply Jte directly later
						mJacScale.de_dsn=-relPose.get_rotation()*  feat_n->getLocalCoordinates();
						list_jacobian_scale.push_back(mJacScale);
					}

				}
				
			}
		}
		std::cout<<"residue/	nb_residue = "<<	residue/	nb_residue<<std::endl;			
		
		
		//accumulate jacobians
		int nb_params=7*nbOptimKf;//scale + translation + rotation
		VectorXf Jte(nb_params);Jte.setZero();
		MatrixXf H(nb_params,nb_params);H.setZero();
				
		//std::cout<<"Update Matrices using jacobianScale"<<std::endl;
		for(int j=0;j<list_jacobian_scale.size();j++)
		{
			kf_loc_jacobian2 &fJacobian=list_jacobian_scale[j];
			
			Jte[7*fJacobian.opt_idn]+=fJacobian.weight * fJacobian.de_dsn.transpose()*fJacobian.pos_error;
			Jte.segment(7*fJacobian.opt_idn+1,6)+=fJacobian.weight * fJacobian.de_dpn.transpose()*fJacobian.pos_error;
			
			H(7*fJacobian.opt_idn,7*fJacobian.opt_idn)+=fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsn;				
			H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpn;				
		}
		//keep information matrix wrt relative pose
		MatrixXf InformationMatrixRelPose=H.block(1,1,6,6);
		
		Jte.segment(1,6)=dHErr_dpn.transpose()*Jte.segment(1,6);
		H.block(1,1,6,6)=dHErr_dpn.transpose()*H.block(1,1,6,6)*dHErr_dpn.transpose();

		Eigen::FullPivLU<MatrixXf> lu(H);
		MatrixXf invH=lu.inverse();
		
		if(isnan(invH(0,0)))
		{
			std::cout<<"invHessian is nan"<<std::endl;
			break;
		}
		else
		{
			VectorXf Dp(nb_params);
			Dp=-gain*(invH*Jte);
			
			if(verb_BA)std::cout<<"\tapply update"<<std::endl;
			float ds=Dp[0];
			//scales_kfs[j]=(1.+ds)*scales_kfs[j];
			scales_kf2=ds+scales_kf2;
			VectorXf dp(Dp.segment(1,6));
			pose_kf2=HomogeneousMatrix22(dp)* pose_kf2;
			
		}
		
		optRelScale=scales_kf2;
		optRelPose=pose_c*pose_n.inverse();
		infoScale=H(0,0);
		infoPose=InformationMatrixRelPose;		
	}
	
	//for testing
	/*for(int f=0;f<KFn.getNbLocalBestFeatures();f++)
		KFn.getPtLocalBestFeatures(f)->depthInRef*=scales_kf2;
	
	HomogeneousMatrix relBestPose= KFn.getBestRelPose();
	relBestPose.set_translation(relBestPose.get_translation()*scales_kf2);
	KFn.setBestRelPose(relBestPose);*/
	

	
}