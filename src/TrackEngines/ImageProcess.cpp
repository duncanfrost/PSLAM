#include "ImageProcess.h"

bool isInImage(Vector2i imageSize,Vector2f p)
{
	return (p[0]>=0 && p[1]>=0 && p[0]<imageSize[0] && p[1]<imageSize[1]);
}
bool isInImageMargin(Vector2i imageSize,int margin,Vector2f p)
{
	return (p[0]>=margin && p[1]>=margin && p[0]<imageSize[0]-margin && p[1]<imageSize[1]-margin);
}
bool isInImage(cv::Size imageSize,Vector2f p)
{
	return (p[0]>=0 && p[1]>=0 && p[0]<imageSize.width && p[1]<imageSize.height);
}
bool isInImageMargin(cv::Size imageSize,int margin,Vector2f p)
{
	return (p[0]>=margin && p[1]>=margin && p[0]<imageSize.width-margin && p[1]<imageSize.height-margin);
}
bool isInImage(cv::Mat &img,Vector2f p)
{
	return (p[0]>=0 && p[1]>=0 && p[0]<img.size().width && p[1]<img.size().height);
}
bool isInImageMargin(cv::Mat &img,int margin,Vector2f p)
{
	return (p[0]>=margin && p[1]>=margin && p[0]<img.size().width-margin && p[1]<img.size().height-margin);
}

void imageGradients(cv::Mat &Img,cv::Mat &gradxImg,cv::Mat &gradyImg)
{
	if(Img.type()!=CV_8UC1)
	{
		std::cerr<<"imageGradients: input image has to be Black and white"<<std::endl;
	}
	else
	{
		GetGradX(Img,gradxImg);
		GetGradY(Img,gradyImg);
		//GetGradXnoOptim(Img,gradxImg);
		//GetGradYnoOptim(Img,gradyImg);
	}
}
void GetGradX(cv::Mat &I, cv::Mat& I2)
{
	I2.create(I.size().height, I.size().width, CV_32FC1);
	int height    = I.size().height;
	int width     = I.size().width;
	uchar *data1      = (uchar *)I.data;
	float *data2      = (float *)I2.data;
	
	//BwImage  Iwrap(I);
	BwImageFloat  dIxwrap(I2);
	
	#pragma omp parallel for
	for(int i=3;i<height*width-3;i++)
	{
		data2[i]=(2047.0 *(data1[i+1] - data1[i-1])
			+913.0 *(data1[i+2] - data1[i-2])
			+112.0 *(data1[i+3] - data1[i-3]))/8418.0;
	}
	for (unsigned int i2=0 ; i2 < I.size().height ; i2++)
	{
		for (unsigned int j=0 ; j < 3 ; j++)
			dIxwrap[i2][j]=0;
		//for (unsigned int j=3 ; j < I.size().width-3 ; j++)
		//	dIxwrap[i2][j]=derivXF(Iwrap,i2,j);
		for (unsigned int j=I.size().width-3 ; j < I.size().width ; j++)
			dIxwrap[i2][j]=0;
	}
	
}

void GetGradXnoOptim(cv::Mat &I, cv::Mat& dIx)
{
	dIx.create(I.size().height, I.size().width, CV_32FC1);
	BwImage  Iwrap(I);
	BwImageFloat  dIxwrap(dIx);
	
	for (unsigned int i=0 ; i < I.size().height ; i++)
	{
		for (unsigned int j=0 ; j < 3 ; j++)
			dIxwrap[i][j]=0;
		for (unsigned int j=3 ; j < I.size().width-3 ; j++)
			dIxwrap[i][j]=derivXF(Iwrap,i,j);
		for (unsigned int j=I.size().width-3 ; j < I.size().width ; j++)
			dIxwrap[i][j]=0;
	}
}
inline double derivXF(BwImage & fr, int r, int c)
{
	return (2047.0 *(fr[r][c+1] - fr[r][c-1])
			+913.0 *(fr[r][c+2] - fr[r][c-2])
			+112.0 *(fr[r][c+3] - fr[r][c-3]))/8418.0;
}
void GetGradY(cv::Mat &I, cv::Mat& I2)
{

	I2.create(I.size().height, I.size().width, CV_32FC1);
	int height    = I.size().height;
	int width     = I.size().width;
	uchar *data1      = (uchar *)&I.at<unsigned char>(0,0);
	float *data2      = (float *)&I2.at<float>(0,0);
	
	//BwImage  Iwrap(I);
	BwImageFloat  dIxwrap(I2);
	
	#pragma omp parallel for
	for(int i=3*width;i<height*width-3*width;i++)
	{
		data2[i]=(2047.0 *(data1[i+width] - data1[i-width])
			+913.0 *(data1[i+2*width] - data1[i-2*width])
			+112.0 *(data1[i+3*width] - data1[i-3*width]))/8418.0;
	}
	for (unsigned int j=0 ; j < 3 ; j++)
		for (unsigned int i2=0 ; i2 < I.size().width ; i2++)
			dIxwrap[j][i2]=0;
	for (unsigned int j=I.size().height-3 ; j < I.size().height ; j++)
		for (unsigned int i2=0 ; i2 < I.size().width ; i2++)
			dIxwrap[j][i2]=0;
	
}
void GetGradYnoOptim(cv::Mat &I, cv::Mat& dIy)
{
	dIy.create(I.size().height, I.size().width, CV_32FC1);
	BwImage  Iwrap(I);
	BwImageFloat  dIywrap(dIy);
	for (unsigned int i=0 ; i < 3 ; i++)
		for (unsigned int j=0 ; j < I.size().width ; j++)
			dIywrap[i][j]=0;
	for (unsigned int i=3 ; i < I.size().height-3 ; i++)
		for (unsigned int j=0 ; j < I.size().width ; j++)
			dIywrap[i][j]=derivYF(Iwrap,i,j);
	for (unsigned int i=I.size().height-3 ; i < I.size().height ; i++)
		for (unsigned int j=0 ; j < I.size().width ; j++)
			dIywrap[i][j]=0;
		
}

double derivYF(BwImage & fr, int r, int c)
{
	return (2047.0 *(fr[r+1][c] - fr[r-1][c])
			+913.0 *(fr[r+2][c] - fr[r-2][c])
			+112.0 *(fr[r+3][c] - fr[r-3][c]))/8418.0;
}


void GetGauss(cv::Mat &I, cv::Mat& GI)
{
	GI.create((int)(I.size().height/2), (int)(I.size().width/2), CV_8UC1);
	cv::Mat GIx ;
	GetGaussX(I, GIx);
	GetGaussY(GIx, GI);	
}
void GetGaussX(cv::Mat &I, cv::Mat& dIx)
{
	dIx.create(I.size().height,(int)( I.size().width/2), CV_8UC1);
	dIx.setTo(0);
	BwImage  wI(I);
	BwImage  wgI(dIx);

	for (unsigned int i=0 ; i < I.size().height ; i++)
	{
		for (unsigned int j=3 ; j < I.size().width-2 ; j+=2)
		{
			wgI[i][j/2]=GaussXF(wI,i,j);
			//dIx[i][j/2]=I[i][j];
		}
	}			
		
}
double GaussXF(BwImage & fr, int r, int c)
{
	return (1.0*(fr[r][c+1] + fr[r][c-1])+
			2.0 *(fr[r][c]))/4.0;
}
void GetGaussY(cv::Mat &I, cv::Mat& dIy)
{
	dIy.create((int)( I.size().height/2),I.size().width, CV_8UC1);
	dIy.setTo(0);
	BwImage  wI(I);
	BwImage  wgI(dIy);
	for (unsigned int i=3 ; i < I.size().height-2 ; i+=2)
	{
		for (unsigned int j=0 ; j < I.size().width ; j++)
		{
			wgI[i/2][j]=GaussYF(wI,i,j);
		}
	}			
		
}
double GaussYF(BwImage & fr, int r, int c)
{
	return (1.0*(fr[r+1][c] + fr[r-1][c])+
			2.0 *(fr[r][c]))/4.0;
}
void extractPatch(cv::Mat &Img,Vector2f pixPos,cv::Mat &Patch,int patchSize)
{
	int halfPatchSize;
	halfPatchSize=(patchSize%2==0)?patchSize/2:(patchSize-1)/2;
	
	/*if(pixPos[0]>=halfPatchSize && pixPos[0]<=Img.size().width-1-halfPatchSize && pixPos[1]>=halfPatchSize && pixPos[1]<=Img.size().height-1-halfPatchSize)
	{
		Patch.create(patchSize, patchSize, CV_8UC1);		
		for(int x=0;x<patchSize;x++)
			for(int y=0;y<patchSize;y++)
				//Patch.at<unsigned char>(y,x)=Img.at<unsigned char>(pixPos[1]-halfPatchSize+y,pixPos[0]-halfPatchSize+x);
				Patch.at<unsigned char>(y,x)=Img.at<unsigned char>((int)(pixPos[1]+0.5)-halfPatchSize+y,(int)(pixPos[0]+0.5)-halfPatchSize+x);
	
	}*/
	/*std::cout<<"extract patch "<<pixPos.transpose()<<std::endl;
	std::cout<<"patchSize "<<patchSize<<std::endl;
	std::cout<<"img size "<<Img.size()<<std::endl;*/
		Patch.create(patchSize, patchSize, CV_8UC1);		
		for(int x=0;x<patchSize;x++)
			for(int y=0;y<patchSize;y++)
			{
				short x2=(int)(pixPos[0]+0.5)-halfPatchSize+x;
				short y2=(int)(pixPos[1]+0.5)-halfPatchSize+y;
				//Patch.at<unsigned char>(y,x)=Img.at<unsigned char>(pixPos[1]-halfPatchSize+y,pixPos[0]-halfPatchSize+x);
				if(isInImage(Img,Vector2f(x2,y2)))
				{
				  	//std::cout<<x2<<" , "<<y2<<std::endl;
					Patch.at<unsigned char>(y,x)=Img.at<unsigned char>(y2,x2);
				}
			}
	
	//std::cout<<"extract patch end"<<std::endl;
	
}
double LevelZeroPos(double dLevelPos, int nLevel)
{
  return (dLevelPos + 0.5) * ScaleLevel(nLevel) - 0.5;
}


// 1-D transform from level zero to level N:
double LevelNPos(double dRootPos, int nLevel)
{
  return (dRootPos + 0.5) / ScaleLevel(nLevel) - 0.5;
}
Vector2f ScaleCurrentToZeroLevel(Vector2f p,int _l)
{
	Vector2f p2;
	p2[0]=LevelZeroPos(p[0],_l);
	p2[1]=LevelZeroPos(p[1],_l);
	return p2;
}
Vector2f ScaleZeroToCurrentLevel(Vector2f p,int _l)
{
	Vector2f p2;
	p2[0]=LevelNPos(p[0],_l);
	p2[1]=LevelNPos(p[1],_l);
	return p2;
}
float ScaleLevel(int _l)
{
	if(_l==0)return 1.;
	return pow(2.,_l);
}
float invScaleLevel(int _l)
{
	if(_l==0)return 1.;
	return pow(0.5,_l);
}
/*float ScaleCurrentToZeroLevel(int _l)
{
	if(_l==0)return 1.;
	return pow(2.,_l);
}
float ScaleZeroToCurrentLevel(int _l)
{
	if(_l==0)return 1.;
	return pow(0.5,_l);
}*/
void drawLine(cv::Mat &imrgb,Vector2f p0,Vector2f p1,Vector3f color)
{
	float norm=sqrt((p1-p0).squaredNorm());
	Vector2f lineUnity=(p1-p0)/norm;
	for(int i=0;i<norm;i++)
	{
		Vector2f pdraw=p0+lineUnity*i;
		imrgb.at<cv::Vec3b>((int)(pdraw[1]+0.5),(int)(pdraw[0]+0.5))=cv::Vec3b(color[2],color[1],color[0]);
	}
}
