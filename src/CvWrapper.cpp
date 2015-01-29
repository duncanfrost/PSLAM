//definition of runnable object

#include "CvWrapper.h"


float getColorSubpix(const cv::Mat& img, Vector2f pt)
{
    assert(!img.empty());
    assert(img.channels() == 1);

    int x = (int)pt[0];
    int y = (int)pt[1];

    int x0 = cv::borderInterpolate(x,   img.cols, cv::BORDER_REFLECT_101);
    int x1 = cv::borderInterpolate(x+1, img.cols, cv::BORDER_REFLECT_101);
    int y0 = cv::borderInterpolate(y,   img.rows, cv::BORDER_REFLECT_101);
    int y1 = cv::borderInterpolate(y+1, img.rows, cv::BORDER_REFLECT_101);

    float a = pt[0] - (float)x;
    float c = pt[1] - (float)y;

    //return (uchar)cvRound((img.at<unsigned char>(y0, x0) * (1.f - a) + img.at<unsigned char>(y0, x1) * a) * (1.f - c)
    //                       + (img.at<unsigned char>(y1, x0) * (1.f - a) + img.at<unsigned char>(y1, x1) * a) * c);
    return (((float)img.at<unsigned char>(y0, x0) * (1.f - a) + (float)img.at<unsigned char>(y0, x1) * a) * (1.f - c)
                           + ((float)img.at<unsigned char>(y1, x0) * (1.f - a) + (float)img.at<unsigned char>(y1, x1) * a) * c);

}
float getColorSubpixf(const cv::Mat& img, Vector2f pt)
{
    assert(!img.empty());
    assert(img.channels() == 1);

    int x = (int)pt[0];
    int y = (int)pt[1];

    int x0 = cv::borderInterpolate(x,   img.cols, cv::BORDER_REFLECT_101);
    int x1 = cv::borderInterpolate(x+1, img.cols, cv::BORDER_REFLECT_101);
    int y0 = cv::borderInterpolate(y,   img.rows, cv::BORDER_REFLECT_101);
    int y1 = cv::borderInterpolate(y+1, img.rows, cv::BORDER_REFLECT_101);

    float a = pt[0] - (float)x;
    float c = pt[1] - (float)y;

    return (img.at<float>(y0, x0) * (1.f - a) + img.at<float>(y0, x1) * a) * (1.f - c) + (img.at<float>(y1, x0) * (1.f - a) + img.at<float>(y1, x1) * a) * c;

}
cv::Vec3b getColorSubpixRGB(const cv::Mat& img, cv::Point2f pt)
{
    assert(!img.empty());
    assert(img.channels() == 3);

    int x = (int)pt.x;
    int y = (int)pt.y;

    int x0 = cv::borderInterpolate(x,   img.cols, cv::BORDER_REFLECT_101);
    int x1 = cv::borderInterpolate(x+1, img.cols, cv::BORDER_REFLECT_101);
    int y0 = cv::borderInterpolate(y,   img.rows, cv::BORDER_REFLECT_101);
    int y1 = cv::borderInterpolate(y+1, img.rows, cv::BORDER_REFLECT_101);

    float a = pt.x - (float)x;
    float c = pt.y - (float)y;

    uchar b = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[0] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[0] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[0] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[0] * a) * c);
    uchar g = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[1] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[1] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[1] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[1] * a) * c);
    uchar r = (uchar)cvRound((img.at<cv::Vec3b>(y0, x0)[2] * (1.f - a) + img.at<cv::Vec3b>(y0, x1)[2] * a) * (1.f - c)
                           + (img.at<cv::Vec3b>(y1, x0)[2] * (1.f - a) + img.at<cv::Vec3b>(y1, x1)[2] * a) * c);

    return cv::Vec3b(b, g, r);
}
#include <fstream>

cv::Mat getDisplayableImage(cv::Mat& img)
{
	assert(img.type() == CV_32FC1);
	
	cv::Mat imBW;
	imBW.create(img.size().height, img.size().width, CV_8UC1);

	
	//get min and max values
	float min=img.at<float>(0,0);
	float max=min;
	
	for (unsigned int j=0 ; j < img.size().width ; j++)
		for (unsigned int i=0 ; i < img.size().height ; i++)
		{
			if(min<img.at<float>(i,j))min=img.at<float>(i,j);
			if(max>img.at<float>(i,j))max=img.at<float>(i,j);
		}
	
	std::cout<<"min = "<<min<<std::endl;
	std::cout<<"max = "<<max<<std::endl;
	
	for (unsigned int j=0 ; j < img.size().width ; j++)
		for (unsigned int i=0 ; i < img.size().height ; i++)
		{
			imBW.at<uchar>(i,j)=255.*(img.at<float>(i,j)-min)/(max-min);	
		}
	
	
	return imBW;

}
