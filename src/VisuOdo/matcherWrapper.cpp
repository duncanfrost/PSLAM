
#include "matcherWrapper.h"

using namespace cv;

matcherWrapper::matcherWrapper()
{
	mlvl_viso=lvl_Viso;
	Homography=Matrix3f::Identity();
	useKLTpoints=false;
	useKLTpointsOnly=false;
	
	MAX_COUNT_KLT=500;
	dist_klt_feat=10;
}

void matcherWrapper::InitWithRef(cv::Mat *imgRef)
{
	uint8_t* img_data  = imgRef[mlvl_viso].data;
	
	// image dimensions
	int32_t dims[3];
	dims[0] = imgRef[mlvl_viso].cols; // width
	dims[1] = imgRef[mlvl_viso].rows; // height
	dims[2] = dims[0]; // image width, BW so it equals to width
		
	if(!useKLTpointsOnly || !useKLTpoints)
		matcher.pushBack(img_data,dims,true);	
	else
	{
		Size subPixWinSize(10,10);
		TermCriteria termcrit(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 20, 0.03);
				
		vector<Point2f> pointsRef;
		goodFeaturesToTrack(imgRef[0], pointsRef, MAX_COUNT_KLT, 0.01, dist_klt_feat, Mat(), 3, 0, 0.04);//get good LK corners
		int margin=5+1;//viso margin
		for(int i=0;i<pointsRef.size();i++)
			if(pointsRef[i].x< ScaleLevel(mlvl_viso)*margin || pointsRef[i].x >imgRef[0].size().width-ScaleLevel(mlvl_viso)*margin 
			    || pointsRef[i].x< ScaleLevel(mlvl_viso)*margin || pointsRef[i].x >imgRef[0].size().width-ScaleLevel(mlvl_viso)*margin)
			{
				pointsRef.erase(pointsRef.begin()+i);
				i--;
			}	
		cornerSubPix(imgRef[0], pointsRef, subPixWinSize, Size(-1,-1), termcrit);//refine position of corners
		
		std::vector<p_feat> kltPointsRef;
		for(int i=0;i<pointsRef.size();i++)
		{
			//p_feat newFeat(pointsRef[i].x/ScaleLevel(mlvl_viso),pointsRef[i].y/ScaleLevel(mlvl_viso),0,4);
			p_feat newFeat(LevelNPos(pointsRef[i].x,mlvl_viso),LevelNPos(pointsRef[i].y,mlvl_viso),0,4);

			kltPointsRef.push_back(newFeat);
		}
		
		matcher.pushBackRef(img_data,dims,kltPointsRef);	
	}
}

void matcherWrapper::match(cv::Mat *_img_p)
{
  	Matrix3f Homography_viso=Homography;
	float div_lvl=ScaleLevel(mlvl_viso);
	for(int i=0;i<2;i++)Homography_viso(i,2)=Homography(i,2)/div_lvl;
	for(int i=0;i<2;i++)Homography_viso(2,i)=Homography(2,i)*div_lvl;
	
	cv::Mat img_2_warped = warpImageInv(_img_p[mlvl_viso],Homography_viso);
 
	//init viso structure
	uint8_t* img_data  = img_2_warped.data;
	int32_t dims[3];
	dims[0] = img_2_warped.cols; // width
	dims[1] = img_2_warped.rows; // height
	dims[2] = dims[0]; // image width, BW so it equals to width
	
	//isntead: Use only warped features on current image	
	std::vector<p_feat> feature_current=matcher.getFeatures(_img_p[mlvl_viso].data,dims);
	//add klt ones
	if(useKLTpointsOnly || useKLTpoints)
	{
		Size subPixWinSize(10,10);
		TermCriteria termcrit(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 20, 0.03);
		
		vector<Point2f> pointsCurr;
	  	goodFeaturesToTrack(_img_p[0], pointsCurr, MAX_COUNT_KLT, 0.01, dist_klt_feat, Mat(), 3, 0, 0.04);//get good LK corners
		//remove features out of viso margin
		int margin=5+1;//viso margin
		for(int i=0;i<pointsCurr.size();i++)
			if(pointsCurr[i].x< ScaleLevel(mlvl_viso)*margin || pointsCurr[i].x >_img_p[0].size().width-ScaleLevel(mlvl_viso)*margin 
			    || pointsCurr[i].x< ScaleLevel(mlvl_viso)*margin || pointsCurr[i].x >_img_p[0].size().width-ScaleLevel(mlvl_viso)*margin)
			{
				pointsCurr.erase(pointsCurr.begin()+i);
				i--;
			}
			  
		cornerSubPix(_img_p[0], pointsCurr, subPixWinSize, Size(-1,-1), termcrit);//refine position of corners

		for(int i=0;i<pointsCurr.size();i++)
		{
			//p_feat newFeat(pointsCurr[i].x/ScaleLevel(mlvl_viso),pointsCurr[i].y/ScaleLevel(mlvl_viso),0,4);
			p_feat newFeat(LevelNPos(pointsCurr[i].x,mlvl_viso),LevelNPos(pointsCurr[i].y,mlvl_viso),0,4);
			feature_current.push_back(newFeat);
		}
	}
	
	//need to apply inverse homography on these features
	Eigen::FullPivLU<MatrixXf> lu(Homography_viso);		
	Matrix3f _Hinv=lu.inverse();
	
	std::vector<p_feat> feature_on_warped;
	for(int i=0;i<feature_current.size();i++)
	{
		p_feat warpedFeat=feature_current[i];
		float x=warpedFeat.u1c;
		float y=warpedFeat.v1c;
		
		
		double denom=_Hinv(2,0)*x+_Hinv(2,1)*y+_Hinv(2,2);
		warpedFeat.u1c=(_Hinv(0,0)*x+_Hinv(0,1)*y+_Hinv(0,2))/denom;
		warpedFeat.v1c=(_Hinv(1,0)*x+_Hinv(1,1)*y+_Hinv(1,2))/denom;
		feature_on_warped.push_back(warpedFeat);
	}
	
	
	//now pass that to matcher to extract descriptors on those (ones that fall into ref when warped onto it (that s why we ll need LUT to get back to id in all features and not only those))
	matcher.pushBack(img_data,dims,feature_on_warped);

	//do matching using these priors
	matcher.matchFeatures();  
	if(useKLTpointsOnly)
		matches=matcher.getMatchesFromClass(4);
	else
		matches=matcher.getMatches();
	
	matcher.bucketFeatures(5,50,50); 
	if(useKLTpointsOnly)
		matchesBucket=matcher.getMatchesFromClass(4);
	else
		matchesBucket=matcher.getMatches();

	matcher.toWarpedId(matches);
	matcher.toWarpedId(matchesBucket);
	
	//current matches got there pos u1c and v1c defined on img_2_warped
	//want them on img current => warp back
	for(int i=0;i<matches.size();i++)
	{
		float x=matches[i].u1c;
		float y=matches[i].v1c;
		
		double denom=Homography_viso(2,0)*x+Homography_viso(2,1)*y+Homography_viso(2,2);
		matches[i].u1c=(Homography_viso(0,0)*x+Homography_viso(0,1)*y+Homography_viso(0,2))/denom;
		matches[i].v1c=(Homography_viso(1,0)*x+Homography_viso(1,1)*y+Homography_viso(1,2))/denom;
	}
	
	if(matches.size()>10)
	{
		//Homography_viso=HomographyFromMatchesRANSAC(matches,100);
		HomographyFromMatchesRobustLM(matches,Homography_viso);
	}
	else
		Homography_viso=Matrix3f::Identity();
	

	//too lvl 0
	if(mlvl_viso!=0)
	{
		for(int i=0;i<matches.size();i++)
		{
			matches[i].u1p=LevelZeroPos(matches[i].u1p,mlvl_viso);
			matches[i].v1p=LevelZeroPos(matches[i].v1p,mlvl_viso);
			matches[i].u1c=LevelZeroPos(matches[i].u1c,mlvl_viso);
			matches[i].v1c=LevelZeroPos(matches[i].v1c,mlvl_viso);

		}
	}
	
	
	//get lvl0 warped
	Homography=Homography_viso;
	for(int i=0;i<2;i++)Homography(i,2)=Homography_viso(i,2)*div_lvl;
	for(int i=0;i<2;i++)Homography(2,i)=Homography_viso(2,i)/div_lvl;
}