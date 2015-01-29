//relative map point class

#pragma once

#include "../AmoDefines.h"
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <vector>


using namespace Eigen;


class MapPoint
{
public:
	MapPoint(){used=true;idKFmeas.reserve(10);i1p.reserve(10);weight=0;};
	
	//position
	Vector3f getPosition(){return position;};
	void updatePosition(Vector3f _p){position=_p;};
	
	//Id
	int getId(){return id;};
	void setId(int newId){id=newId;};
	
	//is point good or bad, bad happens when it is outlier in BA => remove it soon
	bool isUsed(){return used;};
	void setAsBad(){used=false;};
	
	//use weight
	float getWeight(){return weight;};
	void setWeight(float _w){weight=_w;};
	
#ifdef SAVE_POINT_COLOR
	//unsigned char getGrayVal(){return grayVal;};
	//void setGrayVal(unsigned char _c){grayVal=_c;};
	unsigned char getCol(int i){return col[i];};
	void setCol(unsigned char _c[3]){for(int i=0;i<3;i++)col[i]=_c[i];};
#endif
	
	//outlier number
	void consideredInlier(){inlierCount++;};
	void consideredOutlier(){outlierCount++;if(outlierCount>20 && outlierCount>inlierCount)used=false;};
		
	//viewing information: which keyframe sees it, how many
	int nbViews(){return idKFmeas.size();};
	void addView(int idkf,int _i1p){idKFmeas.push_back(idkf);i1p.push_back(_i1p);};
	
	int getView(int i){return idKFmeas[i];};
	int getI1p(int i){return i1p[i];};
	bool doesKFviewMe(int kfId){return !(std::find(idKFmeas.begin(), idKFmeas.end(), kfId)==idKFmeas.end());};
	void removeView(int idkf);
	void removeViewId(int idview);
	void removeViews(){idKFmeas.clear();};
	
	//io functions
	void saveToStream(std::ofstream &fout);
	void loadFromStream(std::ifstream &fout);	
protected:
	//is point any good
	bool used;
	//relative position with Keyframe
	Vector3f position;
	//index in keyframe
	int id;
	
	//tracking information
	short outlierCount;//number of time considered outlier by robust optim
	short inlierCount;//number of time considered outlier by robust optim
	
	float weight;//defined using sum of featureScore of each view
	
	//kf it is measured in:
	std::vector<int> idKFmeas;
	//Viso id of feature corresponding to view => makes it faster to find local feature corresponding to view
	std::vector<int> i1p;
	
#ifdef SAVE_POINT_COLOR
	//unsigned char grayVal;
	unsigned char col[3];
#endif
};


