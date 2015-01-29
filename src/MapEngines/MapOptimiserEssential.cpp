#include "MapOptimiserEssential.h"

MapOptimiserEssential::MapOptimiserEssential(obMap *_myMap)
{
	myMap=_myMap;
	gain=DEFAULT_GAIN_GN;
}


struct jacobianEssential
{
	int opt_id;//index of kf in optimisation list
	float error;//projection error
	MatrixXf de_dp;//jacobian wrt kf position
	float weight;
};

struct jacobianScale
{
	int opt_id;//index of kf in optimisation list
	Vector3f error;//projection error
	Vector3f de_ds;//jacobian wrt kf scale
	MatrixXf de_dt;//jacobian wrt kf translation
	float weight;
};

struct jacobianEssentialcn
{
	int opt_idc;//index of kf in optimisation list
	int opt_idn;//index of kf in optimisation list
	float error;//projection error
	float weight;
	MatrixXf de_dpc;//jacobian wrt kf position
	MatrixXf de_dpn;//jacobian wrt kf position
};

struct jacobianScalecn
{
	int opt_idc;//index of kf in optimisation list
	int opt_idn;//index of kf in optimisation list
	Vector3f error;//projection error
	Vector3f de_dsc;//jacobian wrt kf scale
	MatrixXf de_dtc;//jacobian wrt kf translation
	Vector3f de_dsn;//jacobian wrt kf scale
	MatrixXf de_dtn;//jacobian wrt kf translation
	float weight;
};

struct jacobianPoseRelcn
{
	int opt_idc;//index of kf in optimisation list
	int opt_idn;//index of kf in optimisation list
	float weight;//index of kf in optimisation list
	VectorXf error;//projection error
	MatrixXf InfoMatrix;//jacobian wrt kf translation
	MatrixXf de_dpc;//jacobian wrt kf translation
	MatrixXf de_dpn;//jacobian wrt kf translation
};



void MapOptimiserEssential::optimiseInnerWindowRobust(std::vector<int> &_innerWindowKFs,int nb_iter)
{
  	bool verb_BA=true;
	if(verb_BA)coutGreen<<"###############################################"<<endlGreen;
	if(verb_BA)coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	
	//WARNING just for testing
	//std::vector<int> outerWindow;
	//outerWindow.push_back(_innerWindowKFs[0]-1);
	
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
	if(verb_BA)
	{
		std::cout<<"innerWindowKFs = "<<std::endl;
		for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
		std::cout<<std::endl;
		std::cout<<"outerWindow = "<<std::endl;
		for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
		std::cout<<std::endl<<std::endl;
	}
	
	int nbOptimKf=innerWindowKFs.size();
	
	Camera *myCam=myMap->getCamera();
		
	if(nbOptimKf>0)
	{

	
		if(verb_BA)std::cout<<"\ttest innerWindowKFs passed"<<std::endl;
		//we are going to optimise over scale of each KF and relative pose between KF
		int LUTidKFtoInnerWindow[myMap->getNbKeyFrames()];
		for(int i=0;i<myMap->getNbKeyFrames();i++)LUTidKFtoInnerWindow[i]=-1;
		for(int i=0;i<nbOptimKf;i++)LUTidKFtoInnerWindow[innerWindowKFs[i]]=i;
		  
		
		float mLMLambda = 0.0001;//before LevMarConstant
		float mdLambdaFactor = 2.0;
	
		if(verb_BA)std::cout<<"\tvariables allocated and LUT init"<<std::endl;
		for(int iter=0;iter<nb_iter;iter++)
		{
			//if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			//if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			if(verb_BA)std::cout<<"\t###############################################"<<std::endl;
			if(verb_BA)std::cout<<"\tIteration "<<iter<<std::endl;
			//need to put in minimisation all the points linked by features of innerWindow and outerWindow
			//that should be all the points defined in this keyframes
			std::vector<jacobianEssentialcn> list_jacobian_essential;
			std::vector<jacobianScalecn> list_jacobian_scale;
			
			//residueEssential just to check if error goes down
			float nbEssentialError=0;
			float residueEssential=0;
			float nbScaleError=0;
			float residueScale=0;
			int p_disp=0;//only for debugging
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
					if(std::find(FullWindow.begin(), FullWindow.end(), neigbor.neighboring_kf)!=FullWindow.end())
					{
						if(verb_BA)std::cout<<"\t\tidKFn "<<neigbor.neighboring_kf<<" Nb points : "<<KFn.getNbMapPoint()<<std::endl;
						if(verb_BA)std::cout<<"\t\tidKFn "<<neigbor.neighboring_kf<<" Nb Feat : "<<KFn.getNbLocalBestFeatures()<<std::endl;
						if(verb_BA)std::cout<<"\t\tget Essential Matrix"<<std::endl;
					  
						//compute essential matrix
						HomogeneousMatrix22 pose_c=KFc.getPose();
						HomogeneousMatrix22 pose_n=KFn.getPose();
						
						//HomogeneousMatrix relPose=neigbor.relative_poses;
						HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
						Vector3f t_nc=relPose.get_translation();
						Vector3f nt_nc=t_nc/sqrt(t_nc.squaredNorm());
						Matrix3f R_nc=relPose.get_rotation();
						//Xc=relPose * Xn;
						Matrix3f Ecn=GetSkew(nt_nc)*relPose.get_rotation();
						
						//get Tukey factors
						if(verb_BA)std::cout<<"\t\tget TukeyScale"<<std::endl;
						std::vector<float> vdErrorsEssential;
						std::vector<float> vdErrorsScale;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							//if(verb_BA)std::cout<<"\t\tmatches "<<m<<std::endl;
							p_match &mMatch=neigbor.matches[m];
							//measure in current KF KFc
							//if(verb_BA)std::cout<<"\t\tcompute error"<<std::endl;
							Vector2f measure_n=myCam->ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
							//measure in neigbor
							Vector2f measure_c=myCam->ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
							
							Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
							Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
							
							float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
							//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
							vdErrorsEssential.push_back(errorEssential*errorEssential);

							
							//if(verb_BA)std::cout<<"\t\tcheck scale info "<<m<<std::endl;
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							
								
							uptoscaleFeature *feat_c;
							uptoscaleFeature *feat_n;
							Vector3f scaleError;
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
								//if(verb_BA)std::cout<<"\t\tseems good "<<m<<std::endl;
								feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
								feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
								/*float detph1,recAngle;
								int recGood=reconstructionFromRays(feat_n->posRef,feat_c->posRef,relPose,detph1,recAngle);
								if(recGood!=-1)
								{
									Vector3f feat_n_c=relPose*(feat_n->getLocalCoordinates());
									if(feat_n_c[2]>0);
									{
										scaleError=feat_c->getLocalCoordinates()-feat_n_c;
										vdErrorsScale.push_back(scaleError.squaredNorm());
									}
								}*/
								Vector3f feat_n_c=relPose*(feat_n->getLocalCoordinates());
								scaleError=feat_c->getLocalCoordinates()-feat_n_c;
								vdErrorsScale.push_back(scaleError.squaredNorm());

							}
						}
						if(verb_BA)std::cout<<"\t\tvdErrorsEssential.size() "<<vdErrorsEssential.size()<<std::endl;
						if(verb_BA)std::cout<<"\t\vdErrorsScale.size() "<<vdErrorsScale.size()<<std::endl;
						float sigmaTukeyEssential=getSigmaSquared(vdErrorsEssential);
						float sigmaTukeyScale;
						if(vdErrorsScale.size()!=0)
							sigmaTukeyScale=getSigmaSquared(vdErrorsScale);
						else
							sigmaTukeyScale=0;
						
						
						if(verb_BA)std::cout<<"\tnb matches = "<<neigbor.matches.size()<<std::endl;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							p_match &mMatch=neigbor.matches[m];
							//measure in current KF KFc
							Vector2f measure_n=myCam->ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
							//measure in neigbor
							Vector2f measure_c=myCam->ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
							
							Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
							Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
							
							float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
							//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
							//float TukeyEssCoef=squareRootTukey(errorEssential*errorEssential,sigmaTukeyEssential);
							float TukeyEssCoef=1.;
							residueEssential+=TukeyEssCoef*errorEssential*errorEssential;
							nbEssentialError+=TukeyEssCoef;
							
							
														
							jacobianEssentialcn newJac;
							newJac.opt_idc=-1;
							newJac.opt_idn=-1;
							newJac.error=errorEssential;
							newJac.weight=TukeyEssCoef;
							
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							
							//if(verb_BA)std::cout<<"\t\tmeasure_c = "<<Vector2f(mMatch.u1p,mMatch.v1p).transpose()<<"  tmeasure_n = "<<Vector2f(mMatch.u1c,mMatch.v1c).transpose()<<std::endl;
							//if(verb_BA)std::cout<<"\t\tmMatch.i1p = "<<mMatch.i1p<<" => featc "<<idFeat_c<<"  mMatch.i1c = "<<mMatch.i1c<<" => featn "<<idFeat_n<<std::endl;

								
							uptoscaleFeature *feat_c;
							uptoscaleFeature *feat_n;
							Vector3f scaleError;
							/*float TukeyScaleCoef;
							bool validIntersection=true;
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
								feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
								feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
								
								//check if rays are intersecting (if they are not should not be use for rescaling for sure)
								float detph1,recAngle;
								int recGood=reconstructionFromRays(feat_n->posRef,feat_c->posRef,relPose,detph1,recAngle);
								if(recGood==-1)validIntersection=false;//if ray intersect not properly
								
								Vector3f feat_n_c=relPose*(feat_n->getLocalCoordinates());
								if(feat_n_c[2]<0)validIntersection=false;//if depth of transformed point is negative then probably not well positioned yet for scale
								
								
								scaleError=feat_c->getLocalCoordinates()-feat_n_c;
								TukeyScaleCoef=squareRootTukey(scaleError.squaredNorm(),sigmaTukeyScale);
								residueScale+=TukeyScaleCoef*scaleError.squaredNorm();
								nbScaleError+=TukeyScaleCoef;
							}*/
							
							float TukeyScaleCoef;
							bool validIntersection=true;
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
								//if(verb_BA)std::cout<<"\t\tseems good "<<m<<std::endl;
								feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
								feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
								Vector3f feat_n_c=relPose*feat_n->getLocalCoordinates();
								scaleError=feat_c->getLocalCoordinates()-feat_n_c;
								
								TukeyScaleCoef=squareRootTukey(scaleError.squaredNorm(),sigmaTukeyScale);
								residueScale+=TukeyScaleCoef*scaleError.squaredNorm();
								nbScaleError+=TukeyScaleCoef;
							}
							
							jacobianScalecn mJacScale;
							mJacScale.opt_idc=-1;
							mJacScale.opt_idn=-1;
							mJacScale.error=scaleError;
							if(idFeat_c!=-1 && idFeat_n!=-1)
								mJacScale.weight=TukeyScaleCoef;
								//mJacScale.weight=TukeyScaleCoef*feat_c->recAngle*feat_n->recAngle;
							
							//deriv essential wrt normalised translation
							MatrixXf jac_dnt(1,3);jac_dnt=-hmeasure_c.transpose()*GetSkew(R_nc.transpose()*hmeasure_n);
							Matrix3f dnt_dt=differentiateNormalisedVectorByItself(t_nc);
							MatrixXf jac_dt(1,3);jac_dt=jac_dnt*dnt_dt;
							
							//deriv essential wrt rotation
							MatrixXf A(1,3);A=hmeasure_c.transpose()*GetSkew(nt_nc);
							MatrixXf B(3,1);B=R_nc*hmeasure_n;
							MatrixXf dA_dr(3,3);dA_dr=-GetSkew(hmeasure_c)*GetSkew(nt_nc);
							MatrixXf dB_dr(3,3);dB_dr=-GetSkew(R_nc*hmeasure_n);
							MatrixXf jac_dr(1,3);jac_dr=A*dB_dr+B.transpose()*dA_dr;
							
							//full jacobian wrt se3 transfo nc
							MatrixXf jac_dp(1,6);
							jac_dp.block(0,0,1,3)=jac_dt;
							jac_dp.block(0,3,1,3)=jac_dr;
							//has now to change it to jacobian wrt n and c => need dp_nc/dp_c and dp_nc/dp_n
							
							MatrixXf M1(6,6);
							VectorXf logRelPose=relPose.get_p();
							Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logRelPose[j];
							Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logRelPose[j+3];
							M1.block(0,0,3,3)=-GetSkew(Dw);
							M1.block(3,3,3,3)=-GetSkew(Dw);
							M1.block(0,3,3,3)=-GetSkew(Dt);
							M1.block(3,0,3,3).setZero();
							

							//check if corresponding frames have to be optimised
							if(id_opt_c!=-1)
							{
								//std::cout<<"id_opt_c!=-1"<<std::endl;
								if(TukeyEssCoef>0)
								{
									MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6)+0.5*M1);
									
									newJac.opt_idc=id_opt_c;
									newJac.de_dpc=jac_dp*dHErr_dpc;
								}
								
								if(idFeat_c!=-1 && idFeat_n!=-1 && mJacScale.weight>0)
								{
									mJacScale.opt_idc=id_opt_c;
									//mJacScale.de_dsc=feat_c->getLocalCoordinates();
									mJacScale.de_dsc=relPose*feat_n->getLocalCoordinates();
									mJacScale.de_dtc=-Matrix3f::Identity();
								}
								
								
							}
							if(id_opt_n!=-1)
							{
								if(TukeyEssCoef>0)
								{
									MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1);
									newJac.opt_idn=id_opt_n;
									newJac.de_dpn=jac_dp*dHErr_dpn;
								}
								
								if(idFeat_c!=-1 && idFeat_n!=-1 && TukeyScaleCoef>0)
								{
									mJacScale.opt_idn=id_opt_n;
									//mJacScale.de_dsn=-relPose.get_rotation()*  feat_n->getLocalCoordinates();
									mJacScale.de_dsn=-relPose.get_rotation()*  feat_n->getLocalCoordinates();
									mJacScale.de_dtn=relPose.get_rotation();
								}
							}
							if(id_opt_c!=-1 || id_opt_n!=-1)
							{
								if(newJac.weight>0)
									list_jacobian_essential.push_back(newJac);
								if(idFeat_c!=-1 && idFeat_n!=-1 && validIntersection)
									if(mJacScale.weight>0)
										list_jacobian_scale.push_back(mJacScale);
							}
						}
					}
				}
								
			}
			if(verb_BA)std::cout<<"nbEssentialError = "<<nbEssentialError<<std::endl;
			if(verb_BA)std::cout<<"residueEssential = "<<residueEssential<<std::endl;
			if(verb_BA)std::cout<<"nbScaleError = "<<nbScaleError<<std::endl;
			if(verb_BA)std::cout<<"residueScale = "<<residueScale<<std::endl;
			
			//accumulate jacobians
			int nb_params=7*nbOptimKf;//scale + translation + rotation
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			//std::cout<<"Update Matrices using jacobianEssential"<<std::endl;
			for(int j=0;j<list_jacobian_essential.size();j++)
			{
				jacobianEssentialcn &fJacobian=list_jacobian_essential[j];
				if(fJacobian.opt_idc!=-1)
					Jte.segment(7*fJacobian.opt_idc+1,6)+=fJacobian.weight * fJacobian.de_dpc.transpose()*fJacobian.error;
				if(fJacobian.opt_idn!=-1)
					Jte.segment(7*fJacobian.opt_idn+1,6)+=fJacobian.weight * fJacobian.de_dpn.transpose()*fJacobian.error;
				
				//Hessian:
				if(fJacobian.opt_idc!=-1)
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpc;				
				if(fJacobian.opt_idn!=-1)
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpn;
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1) 
				{
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpn;	
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpc;	
				}
			}	
			
			int NbScaleUpdatePerKf[nbOptimKf];
			for(int i=0;i<nbOptimKf;i++)NbScaleUpdatePerKf[i]=0;
			
			//float lambdaScale=0.001;//make scale less important
			float lambdaScale=1.;//make scale less important
			//std::cout<<"Update Matrices using jacobianScale"<<std::endl;
			for(int j=0;j<list_jacobian_scale.size();j++)
			{
				jacobianScalecn &fJacobian=list_jacobian_scale[j];
				
				if(fJacobian.opt_idc!=-1)
				{
					NbScaleUpdatePerKf[fJacobian.opt_idc]++;
					Jte[7*fJacobian.opt_idc]+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose()*fJacobian.error;
					Jte.segment(7*fJacobian.opt_idc+1,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose()*fJacobian.error;
				}
				if(fJacobian.opt_idn!=-1)
				{
					NbScaleUpdatePerKf[fJacobian.opt_idn]++;
					Jte[7*fJacobian.opt_idn]+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose()*fJacobian.error;
					Jte.segment(7*fJacobian.opt_idn+1,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose()*fJacobian.error;
				}

				
				//Hessian:
				if(fJacobian.opt_idc!=-1)
				{
					H(7*fJacobian.opt_idc,7*fJacobian.opt_idc)+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsc;				
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose() * fJacobian.de_dtc;				
				}
				if(fJacobian.opt_idn!=-1)
				{
					H(7*fJacobian.opt_idn,7*fJacobian.opt_idn)+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsn;				
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose() * fJacobian.de_dtn;				
				}
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1)
				{
					H(7*fJacobian.opt_idc,7*fJacobian.opt_idn)+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsn;				
					H(7*fJacobian.opt_idn,7*fJacobian.opt_idc)+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsc;				
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose() * fJacobian.de_dtc;				
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose() * fJacobian.de_dtn;				
				}
			}
			//check if had enough scale matches to have well conditioned translation Hessian
			//if not have to recondition Hessian cause would have only 2D info from differentiation Essential
			for(int i=0;i<nbOptimKf;i++)
			{
				//coutRed<<"NbScaleUpdatePerKf["<<i<<"] = "<<NbScaleUpdatePerKf[i]<<endlRed;
				if(NbScaleUpdatePerKf[i]==0)
					for(int j=0;j<4;j++)
						H(7*i+j,7*i+j)+=1.;
			}
			
			//std::cout<<"Jte = "<<std::endl;
			//std::cout<<Jte<<std::endl;

			//std::cout<<"H = "<<std::endl;
			//std::cout<<H<<std::endl;
			
			for(int i=0;i<nb_params;i++)
				H(i,i)=(1.+mLMLambda)*H(i,i);

			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			if(isnan(invH(0,0)))
			{
				if(verb_BA)std::cout<<"invHessian is nan"<<std::endl;
				break;
			}
			else
			{
			  
				VectorXf Dp(nb_params);
				Dp=-gain*(invH*Jte);
				
				if(verb_BA)std::cout<<"\tapply update"<<std::endl;
				float scales_kfs[nbOptimKf];
				for(int k=0;k<nbOptimKf;k++)scales_kfs[k]=1;
				HomogeneousMatrix22 pose_kfs[nbOptimKf];
				for(int k=0;k<nbOptimKf;k++)pose_kfs[k]=myMap->getKF(innerWindowKFs[k])->getPose();

				
				//test update to kf
				for(int j=0;j<nbOptimKf;j++)
				{
					if(verb_BA)std::cout<<"\tapply update scale"<<std::endl;
					float ds=Dp[7*j];
					//scales_kfs[j]=(1.+ds)*scales_kfs[j];
					scales_kfs[j]=ds+scales_kfs[j];
					
					if(verb_BA)std::cout<<"\tapply update pose"<<std::endl;
					VectorXf dp(Dp.segment(7*j+1,6));
					pose_kfs[j]=HomogeneousMatrix22(dp)* pose_kfs[j];
				}
				
				
				//recompute residue:
				float nbEssentialErrorAfter=0;
				float residueEssentialAfter=0;
				float nbScaleErrorAfter=0;
				float residueScaleAfter=0;
				for(int k=0;k<FullWindow.size();k++)
				{
					int &idKF=FullWindow[k];
					KeyFrame &KFc=*myMap->getKF(idKF);
					int id_opt_c=LUTidKFtoInnerWindow[idKF];
					
					//get neigbors
					for(int n=0;n<KFc.getNbNeigbours();n++)
					{
						NeigbourKFNew &neigbor=*KFc.getPtNeigbour(n);
						KeyFrame &KFn=*myMap->getKF(neigbor.neighboring_kf);
						int id_opt_n=LUTidKFtoInnerWindow[neigbor.neighboring_kf];
						if(neigbor.neighboring_kf>idKF)//just consider one way, other way should be same error
						{
						  
						  
							//compute essential matrix
							HomogeneousMatrix22 pose_c;
							if(id_opt_c==-1)pose_c=KFc.getPose();
							else pose_c=pose_kfs[id_opt_c];
							
							HomogeneousMatrix22 pose_n;
							if(id_opt_n==-1)pose_n=KFn.getPose();
							else pose_n=pose_kfs[id_opt_n];
							
							//HomogeneousMatrix relPose=neigbor.relative_poses;
							HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
							Vector3f t_nc=relPose.get_translation();
							Vector3f nt_nc=t_nc/sqrt(t_nc.squaredNorm());
							Matrix3f R_nc=relPose.get_rotation();
							//Xc=relPose * Xn;
							//Matrix3f Ecn=GetSkew(t_nc)*relPose.get_rotation();
							Matrix3f Ecn=GetSkew(nt_nc)*relPose.get_rotation();
							
													//get Tukey factors
							std::vector<float> vdErrorsEssential;
							std::vector<float> vdErrorsScale;
							for(int m=0;m<neigbor.matches.size();m++)
							{
								p_match &mMatch=neigbor.matches[m];
								//measure in current KF KFc
								Vector2f measure_n=myCam->ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
								//measure in neigbor
								Vector2f measure_c=myCam->ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
								
								Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
								Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
								
								float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
								//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
								vdErrorsEssential.push_back(errorEssential*errorEssential);

								
								//if matches are linked to local features then give info on scale too
								int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
								int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
								
									
								uptoscaleFeature *feat_c;
								uptoscaleFeature *feat_n;
								Vector3f scaleError;
								if(idFeat_c!=-1 && idFeat_n!=-1)
								{
									feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
									feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
									
									scaleError=feat_c->getLocalCoordinates()-relPose*(feat_n->getLocalCoordinates());
									vdErrorsScale.push_back(scaleError.squaredNorm());
								}
							}
							float sigmaTukeyEssential=getSigmaSquared(vdErrorsEssential);
							float sigmaTukeyScale;
							if(vdErrorsScale.size()!=0)
								sigmaTukeyScale=getSigmaSquared(vdErrorsScale);
							else
								sigmaTukeyScale=0;
							
							for(int m=0;m<neigbor.matches.size();m++)
							{
								p_match &mMatch=neigbor.matches[m];
								//measure in current KF KFc
								Vector2f measure_n=myCam->ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
								//measure in neigbor
								Vector2f measure_c=myCam->ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
								
								Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
								Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
								
								float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
								//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
								//float TukeyEssCoef=squareRootTukey(errorEssential*errorEssential,sigmaTukeyEssential);
								float TukeyEssCoef=1.;
								residueEssentialAfter+=TukeyEssCoef*errorEssential*errorEssential;
								nbEssentialErrorAfter+=TukeyEssCoef;
								
								//if matches are linked to local features then give info on scale too
								int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
								int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
								
								//if(verb_BA)std::cout<<"\t\tmeasure_c = "<<Vector2f(mMatch.u1p,mMatch.v1p).transpose()<<"  tmeasure_n = "<<Vector2f(mMatch.u1c,mMatch.v1c).transpose()<<std::endl;
								//if(verb_BA)std::cout<<"\t\tmMatch.i1p = "<<mMatch.i1p<<" => featc "<<idFeat_c<<"  mMatch.i1c = "<<mMatch.i1c<<" => featn "<<idFeat_n<<std::endl;
								
								//get current estimated scales
								float scale_c;
								if(id_opt_c==-1)scale_c=1;
								else scale_c=scales_kfs[id_opt_c];
									
								float scale_n;
								if(id_opt_n==-1)scale_n=1;
								else scale_n=scales_kfs[id_opt_n];
									
								uptoscaleFeature *feat_c;
								uptoscaleFeature *feat_n;
								Vector3f scaleError;
								if(idFeat_c!=-1 && idFeat_n!=-1)
								{
									feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
									feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
									
									//scaleError=scale_c*feat_c->getLocalCoordinates()-relPose*(scale_n*feat_n->getLocalCoordinates());
									scaleError=feat_c->getLocalCoordinates()-(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
									float TukeyEssCoef=squareRootTukey(scaleError.squaredNorm(),sigmaTukeyScale);
									residueScaleAfter+=TukeyEssCoef*scaleError.squaredNorm();
									nbScaleErrorAfter+=TukeyEssCoef;
								}
							}
						}
					}
									
				}
				if(verb_BA)std::cout<<"nbEssentialErrorAfter = "<<nbEssentialErrorAfter<<std::endl;
				if(verb_BA)std::cout<<"residueEssentialAfter = "<<residueEssentialAfter<<std::endl;
				if(verb_BA)std::cout<<"nbScaleErrorAfter = "<<nbScaleErrorAfter<<std::endl;
				if(verb_BA)std::cout<<"residueScaleAfter = "<<residueScaleAfter<<std::endl;
				
				//if(residueEssentialAfter/nbEssentialErrorAfter<residueEssential/nbEssentialError && residueScaleAfter/nbScaleErrorAfter<residueScale/nbScaleError)
				if(residueEssentialAfter/nbEssentialErrorAfter<residueEssential/nbEssentialError)
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
						
						//VectorXf dp(Dp.segment(7*k+1,6));
						//if(verb_BA)std::cout<<"\t\tpose up "<<HomogeneousMatrix(dp)<<std::endl;
						//if(verb_BA)std::cout<<"\t\tscale up "<<scales_kfs[k]<<std::endl;
						KFc.setPose(pose_kfs[k]);
						for(int f=0;f<KFc.getNbLocalBestFeatures();f++)
							KFc.getPtLocalBestFeatures(f)->depthInRef*=scales_kfs[k];
						
						HomogeneousMatrix22 relPose= KFc.getRelativePose();
						relPose.set_translation(relPose.get_translation()*scales_kfs[k]);
						KFc.setRelativePose(relPose);
						
						HomogeneousMatrix22 relBestPose= KFc.getBestRelPose();
						relBestPose.set_translation(relBestPose.get_translation()*scales_kfs[k]);
						KFc.setBestRelPose(relBestPose);
					}
					
					
				}
				//if not do not confirm change but change lanbda
				else
				{
					mLMLambda = mLMLambda * mdLambdaFactor;
					mdLambdaFactor = mdLambdaFactor * 2;
				}
	
			}
			
			
			
			
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




void MapOptimiserEssential::optimiseInnerWindow(std::vector<int> &_innerWindowKFs,int nb_iter)
{
	coutGreen<<"###############################################"<<endlGreen;
	coutGreen<<"########### optimiseInnerWindow "<<endlGreen;
	innerWindowKFs=_innerWindowKFs;
 	//std::cout<<"innerWindowKFs = "<<std::endl;
	//for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  "<<std::endl;

	//check if we have some exterior window:
	std::vector<int> outerWindow=myMap->getDirectNeigbors(_innerWindowKFs);
	
	//WARNING just for testing
	//std::vector<int> outerWindow;
	//outerWindow.push_back(_innerWindowKFs[0]-1);
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
	std::cout<<"innerWindowKFs = "<<std::endl;
	for(int i=0;i<innerWindowKFs.size();i++)std::cout<<innerWindowKFs[i]<<"  ";
	std::cout<<std::endl;
	std::cout<<"outerWindow = "<<std::endl;
	for(int i=0;i<outerWindow.size();i++)std::cout<<outerWindow[i]<<"  ";
	std::cout<<std::endl<<std::endl;
	
	bool verb_BA=true;
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
			std::vector<jacobianEssentialcn> list_jacobian_essential;
			std::vector<jacobianScalecn> list_jacobian_scale;
			
			//residueEssential just to check if error goes down
			float nbEssentialError=0;
			float residueEssential=0;
			int p_disp=0;//only for debugging
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
					if(std::find(FullWindow.begin(), FullWindow.end(), neigbor.neighboring_kf)!=FullWindow.end())
					{
						//if(verb_BA)std::cout<<"\tuse edge "<<idKF<<" <=> "<<neigbor.neighboring_kf<<std::endl;
						//compute essential matrix error; ie x_i *E_ij *x_j

					  
					  
						//compute essential matrix
						HomogeneousMatrix22 pose_c=KFc.getPose();
						HomogeneousMatrix22 pose_n=KFn.getPose();
						
						//HomogeneousMatrix relPose=neigbor.relative_poses;
						HomogeneousMatrix22 relPose=pose_c*pose_n.inverse();
						Vector3f t_nc=relPose.get_translation();
						Vector3f nt_nc=t_nc/sqrt(t_nc.squaredNorm());
						Matrix3f R_nc=relPose.get_rotation();
						//Xc=relPose * Xn;
						//Matrix3f Ecn=GetSkew(t_nc)*R_nc;
						Matrix3f Ecn=GetSkew(nt_nc)*R_nc;
						
						if(verb_BA)std::cout<<"\tnb matches = "<<neigbor.matches.size()<<std::endl;
						for(int m=0;m<neigbor.matches.size();m++)
						{
							p_match &mMatch=neigbor.matches[m];
							//measure in current KF KFc
							Vector2f measure_n=myCam->ToMeters(Vector2f(mMatch.u1p,mMatch.v1p));
							//measure in neigbor
							Vector2f measure_c=myCam->ToMeters(Vector2f(mMatch.u1c,mMatch.v1c));
							
							Vector3f hmeasure_n(measure_n[0],measure_n[1],1.);
							Vector3f hmeasure_c(measure_c[0],measure_c[1],1.);
							
							//float errorEssential=hmeasure_c.transpose()*GetSkew(t_nc)*R_nc*hmeasure_n;
							float errorEssential=hmeasure_c.transpose()*Ecn*hmeasure_n;
							//std::cout<<"errorEssential = "<<errorEssential<<std::endl;
							residueEssential+=errorEssential*errorEssential;
							nbEssentialError++;
							
							jacobianEssentialcn newJac;
							newJac.opt_idc=-1;
							newJac.opt_idn=-1;
							newJac.error=errorEssential;
							newJac.weight=1.;
							
							//if matches are linked to local features then give info on scale too
							int idFeat_c=KFc.indexCandidateFeatureFromVisoId(mMatch.i1c);
							int idFeat_n=KFn.indexCandidateFeatureFromVisoId(mMatch.i1p);
							
							//if(verb_BA)std::cout<<"\t\tmeasure_c = "<<Vector2f(mMatch.u1p,mMatch.v1p).transpose()<<"  tmeasure_n = "<<Vector2f(mMatch.u1c,mMatch.v1c).transpose()<<std::endl;
							//if(verb_BA)std::cout<<"\t\tmMatch.i1p = "<<mMatch.i1p<<" => featc "<<idFeat_c<<"  mMatch.i1c = "<<mMatch.i1c<<" => featn "<<idFeat_n<<std::endl;
							
							//get current estimated scales
							float scale_c;
							if(id_opt_c==-1)scale_c=1;
							else scale_c=scales_kfs[id_opt_c];
								
							float scale_n;
							if(id_opt_n==-1)scale_n=1;
							else scale_n=scales_kfs[id_opt_n];
								
							uptoscaleFeature *feat_c;
							uptoscaleFeature *feat_n;
							Vector3f scaleError;
							/*bool validIntersection=false;
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
								//if(verb_BA)std::cout<<"\t\tseems good "<<m<<std::endl;
								feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
								feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
								float detph1,recAngle;
								int recGood=reconstructionFromRays(feat_n->posRef,feat_c->posRef,relPose,detph1,recAngle);
								if(recGood!=-1)
								{
									//Vector3f feat_n_c=relPose*(scale_n*feat_n->getLocalCoordinates());
									Vector3f feat_n_c=(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
									if(feat_n_c[2]>0);
									{
										//scaleError=scale_c*feat_c->getLocalCoordinates()-feat_n_c;
										scaleError=feat_c->getLocalCoordinates()-feat_n_c;
										validIntersection=true;
									}
								}
							}*/
							bool validIntersection=true;
							if(idFeat_c!=-1 && idFeat_n!=-1)
							{
								//if(verb_BA)std::cout<<"\t\tseems good "<<m<<std::endl;
								feat_c=KFc.getPtLocalBestFeatures(idFeat_c);
								feat_n=KFn.getPtLocalBestFeatures(idFeat_n);
								Vector3f feat_n_c=(1./scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
								scaleError=feat_c->getLocalCoordinates()-feat_n_c;
							}
							
							jacobianScalecn mJacScale;
							mJacScale.opt_idc=-1;
							mJacScale.opt_idn=-1;
							mJacScale.error=scaleError;
							mJacScale.weight=1.;
							

							//deriv essential wrt normalised translation
							MatrixXf jac_dnt(1,3);jac_dnt=-hmeasure_c.transpose()*GetSkew(R_nc.transpose()*hmeasure_n);
							Matrix3f dnt_dt=differentiateNormalisedVectorByItself(t_nc);
							MatrixXf jac_dt(1,3);jac_dt=jac_dnt*dnt_dt;
							
							//deriv essential wrt rotation
							MatrixXf A(1,3);A=hmeasure_c.transpose()*GetSkew(nt_nc);
							MatrixXf B(3,1);B=R_nc*hmeasure_n;
							MatrixXf dA_dr(3,3);dA_dr=-GetSkew(hmeasure_c)*GetSkew(nt_nc);
							MatrixXf dB_dr(3,3);dB_dr=-GetSkew(R_nc*hmeasure_n);
							MatrixXf jac_dr(1,3);jac_dr=A*dB_dr+B.transpose()*dA_dr;
							
							//full jacobian wrt se3 transfo nc
							MatrixXf jac_dp(1,6);
							jac_dp.block(0,0,1,3)=jac_dt;
							jac_dp.block(0,3,1,3)=jac_dr;
							//has now to change it to jacobian wrt n and c => need dp_nc/dp_c and dp_nc/dp_n
							
							MatrixXf M1(6,6);
							VectorXf logRelPose=relPose.get_p();
							Vector3f Dt;for(int j=0;j<3;j++)Dt[j]=logRelPose[j];
							Vector3f Dw;for(int j=0;j<3;j++)Dw[j]=logRelPose[j+3];
							M1.block(0,0,3,3)=-GetSkew(Dw);
							M1.block(3,3,3,3)=-GetSkew(Dw);
							M1.block(0,3,3,3)=-GetSkew(Dt);
							M1.block(3,0,3,3).setZero();

							
							//check if corresponding frames have to be optimised
							if(id_opt_c!=-1)
							{
								//std::cout<<"id_opt_c!=-1"<<std::endl;
								//MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6));
								MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6)+0.5*M1);
								//MatrixXf dHErr_dpc=(MatrixXf::Identity(6,6)+0.5*M1+M1*M1/12.);

								
								newJac.opt_idc=id_opt_c;
								//newJac.de_dpc=jac_dp;
								newJac.de_dpc=jac_dp*dHErr_dpc;
								
								if(idFeat_c!=-1 && idFeat_n!=-1)
								{
									mJacScale.opt_idc=id_opt_c;
									//mJacScale.de_dsc=scale_c*feat_c->getLocalCoordinates();
									mJacScale.de_dsc=1./(scale_c*scale_c)*(relPose*(scale_n*feat_n->getLocalCoordinates()));
									//std::cout<<"de_dsc = "<<std::endl;
									//std::cout<<mJacScale.de_dsc.transpose()<<std::endl;
									mJacScale.de_dtc=-Matrix3f::Identity();
								}	
							}
							if(id_opt_n!=-1)
							{
								MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1);
								//MatrixXf dHErr_dpn=-(MatrixXf::Identity(6,6)-0.5*M1-M1.transpose()*M1/12.);

								
								newJac.opt_idn=id_opt_n;
								newJac.de_dpn=jac_dp*dHErr_dpn;
								
								if(idFeat_c!=-1 && idFeat_n!=-1)
								{
									mJacScale.opt_idn=id_opt_n;
									//mJacScale.de_dsn=-relPose.get_rotation()*  scale_n*feat_n->getLocalCoordinates();
									mJacScale.de_dsn=-(1./scale_c)*relPose.get_rotation()*  feat_n->getLocalCoordinates();
									mJacScale.de_dtn=relPose.get_rotation();
								}
							}
							if(id_opt_c!=-1 || id_opt_n!=-1)
							{
								list_jacobian_essential.push_back(newJac);
								if(validIntersection && idFeat_c!=-1 && idFeat_n!=-1)
									list_jacobian_scale.push_back(mJacScale);
							}

						}
					}
				}
								
			}
			std::cout<<"nbEssentialError = "<<nbEssentialError<<std::endl;
			std::cout<<"residueEssential = "<<residueEssential<<std::endl;
			
			//accumulate jacobians
			int nb_params=7*nbOptimKf;//scale + translation + rotation
			VectorXf Jte(nb_params);Jte.setZero();
			MatrixXf H(nb_params,nb_params);H.setZero();
			
			std::cout<<"Update Matrices using jacobianEssential"<<std::endl;
			for(int j=0;j<list_jacobian_essential.size();j++)
			{
				jacobianEssentialcn &fJacobian=list_jacobian_essential[j];
				if(fJacobian.opt_idc!=-1)
					Jte.segment(7*fJacobian.opt_idc+1,6)+=fJacobian.weight * fJacobian.de_dpc.transpose()*fJacobian.error;
				if(fJacobian.opt_idn!=-1)
					Jte.segment(7*fJacobian.opt_idn+1,6)+=fJacobian.weight * fJacobian.de_dpn.transpose()*fJacobian.error;
				
				//Hessian:
				if(fJacobian.opt_idc!=-1)
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpc;				
				if(fJacobian.opt_idn!=-1)
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpn;
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1) 
				{
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,6,6)+=fJacobian.weight * fJacobian.de_dpc.transpose() * fJacobian.de_dpn;	
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,6,6)+=fJacobian.weight * fJacobian.de_dpn.transpose() * fJacobian.de_dpc;	
				}
			}		
	
			int NbScaleUpdatePerKf[nbOptimKf];
			for(int i=0;i<nbOptimKf;i++)NbScaleUpdatePerKf[i]=0;
			
			float lambdaScale=0.001;//make scale less important
			std::cout<<"Update Matrices using jacobianScale"<<std::endl;
			for(int j=0;j<list_jacobian_scale.size();j++)
			{
				jacobianScalecn &fJacobian=list_jacobian_scale[j];
				
				if(fJacobian.opt_idc!=-1)
				{
					NbScaleUpdatePerKf[fJacobian.opt_idc]++;
					Jte[7*fJacobian.opt_idc]+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose()*fJacobian.error;
					Jte.segment(7*fJacobian.opt_idc+1,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose()*fJacobian.error;
				}
				if(fJacobian.opt_idn!=-1)
				{
					NbScaleUpdatePerKf[fJacobian.opt_idn]++;
					Jte[7*fJacobian.opt_idn]+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose()*fJacobian.error;
					Jte.segment(7*fJacobian.opt_idn+1,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose()*fJacobian.error;
				}
				
				//Hessian:
				if(fJacobian.opt_idc!=-1)
				{
					H(7*fJacobian.opt_idc,7*fJacobian.opt_idc)+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsc;				
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idc+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose() * fJacobian.de_dtc;				
				}
				if(fJacobian.opt_idn!=-1)
				{
					H(7*fJacobian.opt_idn,7*fJacobian.opt_idn)+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsn;				
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idn+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose() * fJacobian.de_dtn;				
				}
				if(fJacobian.opt_idc!=-1 && fJacobian.opt_idn!=-1)
				{
					H(7*fJacobian.opt_idc,7*fJacobian.opt_idn)+=lambdaScale*fJacobian.weight * fJacobian.de_dsc.transpose() * fJacobian.de_dsn;				
					H(7*fJacobian.opt_idn,7*fJacobian.opt_idc)+=lambdaScale*fJacobian.weight * fJacobian.de_dsn.transpose() * fJacobian.de_dsc;				
					H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose() * fJacobian.de_dtc;				
					H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose() * fJacobian.de_dtn;				
					//H.block(7*fJacobian.opt_idn+1,7*fJacobian.opt_idc+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtc.transpose() * fJacobian.de_dtn;				
					//H.block(7*fJacobian.opt_idc+1,7*fJacobian.opt_idn+1,3,3)+=lambdaScale*fJacobian.weight * fJacobian.de_dtn.transpose() * fJacobian.de_dtc;				
				}

			}
			
			//check if had enough scale matches to have well conditioned translation Hessian
			//if not have to recondition Hessian cause would have only 2D info from differentiation Essential
			for(int i=0;i<nbOptimKf;i++)
			{
				coutRed<<"NbScaleUpdatePerKf["<<i<<"] = "<<NbScaleUpdatePerKf[i]<<endlRed;
				if(NbScaleUpdatePerKf[i]==0)
					H.block(7*i+1,7*i+1,3,3)+=Matrix3f::Identity();
			}


			Eigen::FullPivLU<MatrixXf> lu(H);
			MatrixXf invH=lu.inverse();
			std::cout<<"rank = "<<lu.rank()<<std::endl;
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
			
		}
	}
}
	
}



