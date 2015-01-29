//Keyframe class: the whole map is stored as a collection
//of keyframe. With in each, an image and all the features 
//and descriptors necessary to do any matching or tracking
//and a collection of points and measures. 

//it has been prefered to store the points in their source 
//KF so that if we are interested in very big map then 
//loading a submap of the map would be much faster and 
//simpler. Also it makes it easy to define points with 
//relative poses with respect to KF so that if we perform a
//pose graph optimisation and uptade the position of the 
//KFs then the position of the points would not have to be 
//changed (for a coarse but not too bad solution...)

#pragma once

#include "../AmoDefines.h"
#include <vector>
#include "../Primitives/HomogeneousMatrix.h"
#include "../Primitives/MapPoint.h"
#include "../VisuOdo/matcherWrapper.h"
#include "../TrackEngines/GlobalTransfoEstim.h"
#include "../MapEngines/GeoFunctions.h"
#include "../TrackEngines/ImageProcess.h"
#include "../MapEngines/HomogFinder.h"

#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/imgproc/imgproc.hpp"


class KeyFrame;

struct uptoscaleFeature
{
	//TODO integrate fundamental matric score from which feature has been created, could be good info to have robust rescaling
	Vector2f posRef;//position of feature in meters in KF
	float depthInRef;//estimated depth
	float recAngle;//reconstruction angle	
	float scoreFundamentalOrigin;
	int i1p;//id of corresponding corner
	
	int nb_outlier;
	
#ifdef SAVE_POINT_COLOR
	//unsigned char grayVal;
	unsigned char col[3];
#endif
	
	//link to Map point, ie resulting estimation out of one or
	//more measures
	bool matched;
	//int KForigin;
	KeyFrame *ptKForigin;
	int idPoint;	
	
	uptoscaleFeature()
	{
		matched=false;
		nb_outlier=0;
	};
	
	Vector3f getLocalCoordinates(){return toHomogeneous(posRef)*depthInRef;};
	
		
	//io functions
	void saveToStream(std::ofstream &fout, KeyFrame *kf0);
	void loadFromStream(std::ifstream &fout, KeyFrame *kf0);	
	
};

struct NeigbourKFNew
{
	int neighboring_kf; // id of neigboring KF
	Matrix3f Homography;//homography I_neigbor(pix)=I(H(pix))
	std::vector<p_match> matches;//matches viso between two images (from neigbor to this)
	float edgeScore;//the edge score is equal to sum of (multiplication of fundamental matrix scores of interconnected local features)

	//not very yet....
	float relative_scale;//scale to apply to neigbor to match this
	float Informationscale;
	HomogeneousMatrix relative_poses; //relative pose from neigbor to this using this scale
	MatrixXf InformationMatrixPose;
	
		//io functions
	void saveToStream(std::ofstream &fout)
	{
		fout.write((const char*)&neighboring_kf,sizeof(int));	
		relative_poses.saveToStream(fout);
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)
		      fout.write((const char*)&Homography(i,j),sizeof(float));
		int nb_matches=matches.size();
		fout.write((const char*)&nb_matches,sizeof(int));	
		for(int i=0;i<nb_matches;i++)
			fout.write((const char*)&matches[i],sizeof(p_match));
		fout.write((const char*)&edgeScore,sizeof(float));
		for(int i=0;i<6;i++)for(int j=0;j<6;j++)
		      fout.write((const char*)&InformationMatrixPose(i,j),sizeof(float));
		fout.write((const char*)&relative_scale,sizeof(float));	
		fout.write((const char*)&Informationscale,sizeof(float));	
	  
	}
	void loadFromStream(std::ifstream &fout)
	{
		fout.read((char*)&neighboring_kf,sizeof(int));
		relative_poses.loadFromStream(fout);
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)
		      fout.read((char*)&Homography(i,j),sizeof(float));
		int nb_matches;
		fout.read((char*)&nb_matches,sizeof(int));
		matches.resize(nb_matches);		
		for(int i=0;i<nb_matches;i++)
			fout.read((char*)&matches[i],sizeof(p_match));
		fout.read((char*)&edgeScore,sizeof(float));
		InformationMatrixPose.resize(6,6);
		for(int i=0;i<6;i++)for(int j=0;j<6;j++)
		      fout.read((char*)&InformationMatrixPose(i,j),sizeof(float));
		//std::cout<<"InformationMatrixPose size :"<<InformationMatrixPose.rows()<<" x "<<InformationMatrixPose.cols()<<std::endl;
		//std::cout<<InformationMatrixPose<<std::endl;
		fout.read((char*)&relative_scale,sizeof(float));
		fout.read((char*)&Informationscale,sizeof(float));
	}
};

#define max_view_in_mini_BA 10
struct matchMiniBA
{
	bool matched;
	Vector2f pos_p;//pos in meter in ref
	Vector2f pos_c;//pos in meter in match
	matchMiniBA(){matched=false;};
};

class KeyFrame
{
public:
	KeyFrame();
	void InitMemory();
	KeyFrame(int id,cv::Mat &_img,HomogeneousMatrix _pose);
	void Init(int id,cv::Mat &_img,HomogeneousMatrix _pose);
	~KeyFrame();	
	
	//init matching tools
	void initMatcher();
	void extractORBFeatures();
	
	//matching current image with KF providing coarse rotation transformation between kf and previous image
	int useNewFrame(cv::Mat *_img_c,Camera *_myCamera);
	//estimate 3d features from relative position and set of matches
	std::vector<uptoscaleFeature> getGoodFeaturesFromRayIntersection(std::vector<p_match> &matches,Camera *myCamera,HomogeneousMatrix &KFtoCurrent);	
	//do bundle adjustment to refine algebraic solution
	//for now only use in Testing Virtual Scene (when put considerable image noise algebraic solution is far from refined one)
	void doMiniBA(Camera *myCamera){doMiniBA(best_matches,best_features,myCamera,bestRelPose);};
	void doMiniBA(std::vector<p_match> &matches,std::vector<uptoscaleFeature> &feats, Camera *myCamera,HomogeneousMatrix &poseRel,bool useLM=false);
	//check if depth is consistent in new stereo
	std::vector<uptoscaleFeature> filterWithDepthConsistancy(std::vector<uptoscaleFeature> &GoodFeatures1);
	//get fundamental matrix score and average reconstruction angle of 3d points
	float getFundamentalMatrixScore(std::vector<uptoscaleFeature> &GoodFeatures2,float &_avr_angle);
	//use 3d features
	HomogeneousMatrix RescalePose(HomogeneousMatrix &KFtoCurrentnoscale,std::vector<uptoscaleFeature> &GoodFeatures2,std::vector<uptoscaleFeature> &GoodFeatures1,std::vector<uptoscaleFeature> &best_features);

	//initialisation of local 3d using one neigbor keyframe: match between neigbor and current frame already done
	//as well as 3D relative pose estimation => use inverse of those to init best fundamental and bestFeatures
	void useNeigborForInitLocalStereo(cv::Mat &imgBest,std::vector<p_match> &matches, HomogeneousMatrix &_relPoseNeigb,Camera *_myCamera);
	//void displayMatchc(int i);

	//##################################################################
	//Things that have to do with map:		
	//pose of the keyFrame
	void setPose(HomogeneousMatrix _w_To_cam){w_To_cam.Init(_w_To_cam);};
	HomogeneousMatrix getPose(){return w_To_cam;};
	
	//id 
	void setId(int _id){id=_id;};
	int getId(){return id;};
	
	//called by constructor to init image pyramid
	void makeKF_Image(cv::Mat &_img);	
	cv::Mat *getImg_p(){return img_p;};
	cv::Mat &getImg_p(int _l){return img_p[_l];};
	
	//depth range
	Vector2f getDepthRange(){return depthRange;};
	void setDepthRange(Vector2f _dr){depthRange=_dr;};
	
	//get local best stereo result with scale
	int getNbLocalBestFeatures(){return best_features.size();};
	void addLocalBestFeature(uptoscaleFeature &_c){best_features.push_back(_c);};
	uptoscaleFeature getLocalBestFeatures(int i){return best_features[i];};
	uptoscaleFeature *getPtLocalBestFeatures(int i){return &best_features[i];};
	std::vector<uptoscaleFeature> *getLocalBestFeatures(){return &best_features;};
	void removeLocalBestFeatures(int i){best_features.erase(best_features.begin()+i);};//warning should remove Point view if feature is matched
	void clearLocalBestFeatures(){best_features.clear();};
	int indexCandidateFeatureFromVisoId(int idFeaturep);
	int getNumberOfFeatureMatched();
	HomogeneousMatrix getBestRelPose(){return bestRelPose;};
	void setBestRelPose(HomogeneousMatrix _h){bestRelPose=_h;};
	cv::Mat &getBestImgPair(){return img_best_pair;};
	std::vector<p_match> *getBestLocalMatches(){return &best_matches;};
	void setBestLocalMatches(std::vector<p_match> &_bestM){best_matches=_bestM;};
	
	//use neigborhood
	bool addNeighbour(NeigbourKFNew &_nKF){neigbours.push_back(_nKF);};
	int getNbNeigbours(){return neigbours.size();};
	NeigbourKFNew getNeigbour(int i){return neigbours[i];};//warning quite a large structure since has all matches stored in it, prefer pointer
	NeigbourKFNew *getPtNeigbour(int i){return &neigbours[i];};
	int isNeigbour(int i);
	
	//get map points
	int getNbMapPoint(){return mMapPoints.size();};
	void addMapPoint(MapPoint &_p){mMapPoints.push_back(_p);};
	MapPoint *getPtMapPoint(int i){return &mMapPoints[i];};
	void removePoint(int i){mMapPoints.erase(mMapPoints.begin()+i);};
	int getNbUsedMapPoint(){int res=0;for(int i=0;i<mMapPoints.size();i++)if(mMapPoints[i].isUsed())res++;return res;};
	
	//##################################################################
	//Things that have to do with tracker:	
	//get relative pose estimate with lastely seen image
	HomogeneousMatrix getRelativePose(){return relPose;};
	void setRelativePose(HomogeneousMatrix &_p){relPose=_p;};
	void setRelativePosePrevious(HomogeneousMatrix &_p){relPosePrevious=_p;};

	//estimate relative pose of last matched frame using matches and best stereo
	HomogeneousMatrix computeRelativeCurrentPoseWithLocalFeatures(Camera *_myCamera);		
	HomogeneousMatrix computeRelativeCurrentPoseWithMatchedFeatures(Camera *_myCamera);		
	HomogeneousMatrix computeRelativeCurrentPoseWithAllMatches(Camera *_myCamera);		

	//get last matches
	int getLastNumberMatches(){return matchesCurrent.size();};
	std::vector<p_match> getCurrentMatches(){return matchesCurrent;};
	
	//get overlap with lastly seen frame; just use lastly estimated homography
	void computeOverlap();
	float getOverlapWithLastFrame(){return overlapWithLast;};
	
	//get overlap with lastly seen frame; just use lastly estimated homography
	Matrix3f getHomography(){return Homography;};
	void setHomography(Matrix3f _h){Homography=_h;matcher.setHomography(Homography);};
	
	//get average reconstruction angle with lastly seen frame (used in mapTracker to check if baseline large enough to define new keyframe)
	float getAverageReconstructionAngle(){return average_recAngle;};
	
	//get last estimated fundamental score
	float getLastFundamentalScore(){return last_fundamental_score;}
	
	//see name
	bool haveLocalBestFeaturesChanged(){return localBestFeaturesHaveChanged;};
	void LocalBestFeaturesHaveBeenChecked(){localBestFeaturesHaveChanged=false;};
	
	//void SaveNewMatchesAndPoseForMiniBA(std::vector<p_match> &matchesCurrent,HomogeneousMatrix &relPose,Camera *_myCamera);
	//void DoMiniBA(Camera *_myCamera);
	
	//io functions
	void saveToStream(std::ofstream &fout, KeyFrame *kf0);
	void loadFromStream(std::ifstream &fout, KeyFrame *kf0);	
private:
	//index of frame
	int id;
  
	//bool fixed;
	HomogeneousMatrix w_To_cam;
	
	//current camera rotation
	//Vector3f rotation;
	//Homography between keyframe and lastly seen image; ie I(H(pix))=I_ref(pix)
	Matrix3f Homography;//warning defined in lvl_Viso pyramid level
	//NCC with last image
	float overlapWithLast;
	
	//image corresponding to keyframe in black and white
	cv::Mat img_p[NB_LEVELS];//pyramidal
	cv::Mat img_col;
	
	//for continuous tracking
	//Matcher matcher;
	matcherWrapper matcher;
	std::vector<p_match> matchesCurrent;//matches with last frame given to KF (from ref to last)
	
	//for loop closure
	//list of features in keyframes with corresponding measure and link to mapped 3D point if so
	std::vector<Vector2f> OrbCorners[NB_LEVELS];	
	//descriptors
	cv::Mat OrbDescriptors[NB_LEVELS];
	
	//depth range
	Vector2f depthRange;	
	
	//scale of relative map in KF
	//float scale;
	
	//average reconstruction angle of points and fundamental matrix score with last frame
	float average_recAngle;	
	float last_fundamental_score;
	
	//if had bestFeatures linked to map and have new ones=> have to check with neigbors if new connections are possible
	bool localBestFeaturesHaveChanged; 
	
	//lastely estimated relative pose between current frame and Keyframe, coarse estim from fundamental matrix then rescale, then can be refined
	HomogeneousMatrix relPose;
	//last pose estimated by tracker
	HomogeneousMatrix relPosePrevious;
	
	//info resulting from best stereo
	cv::Mat img_best_pair;
	float best_score_fundamental;
	HomogeneousMatrix bestRelPose;
	std::vector<p_match> best_matches;
	std::vector<uptoscaleFeature> best_features;
	
	//list of connected keyframes by index(sharing enough covisible points in common)
	std::vector<NeigbourKFNew> neigbours;
	
	//list of map points stored in this KF
	std::vector<MapPoint> mMapPoints;
	

};

//express all the variables in frame1//get 3d point from two measures and relative pose
//search for two unknown : alpha and beta so that alpha*v1-beta*v2=c1c2
//with alpha*v1 = c1p (p=intersection) and beta*v2=c2p
//=> c1p+pc2=c1c2
int reconstructionFromRays(Vector2f mes1,Vector2f mes2,HomogeneousMatrix pose1_2,float &depth1, float &recAngle, bool checkIntersect=true);
int reconstructionFromRaysInvDepth(Vector2f mes1,Vector2f mes2,HomogeneousMatrix pose1_2,float &invdepth, float &recAngle);//return 0 if points at infinity

//check if new stereo that we have is better than what we had so far
bool CheckNewStereo(std::vector<p_match> &matches,Camera *_myCamera,HomogeneousMatrix &KFtoCurrent);

//do bundle adjustment to refine algebraic solution
//void doMiniBA(std::vector<p_match> &matches,Camera *myCamera,HomogeneousMatrix &poseRel);
