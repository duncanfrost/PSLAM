//Copyright 2008 Isis Innovation Limited 

#include "ImageWindow.h"

static int nb_Image_windows=0;

namespace ImageWindowNs
{
TextureSet tempTextureImageWindow;
}
  

using namespace ImageWindowNs;

ImageWindow::ImageWindow(std::string _title,amoVideoSource *_vid_source)
{
	title=_title;
	vid_source=_vid_source;
	width=vid_source->getSize()[0];
	height=vid_source->getSize()[1];
	
	pt_img=NULL;
	nb_Image_windows++;
	//vid_source2=_vid_source;
}
ImageWindow::ImageWindow(std::string _title,cv::Mat *_img)
{
	title=_title;
	pt_img=_img;
	width=pt_img->cols;
	height=pt_img->rows;

	vid_source=NULL;
	nb_Image_windows++;
}

void ImageWindow::prepare_draw()
{
	cv::Mat *ptmimFrame;
	if(pt_img==NULL)
		ptmimFrame=vid_source->GetFramePointer();
	else
		ptmimFrame=pt_img;
	LoadTexture(*ptmimFrame,tempTextureImageWindow);

	//use to be in drawImageWindow but was causing a synchronisation problem
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	DrawTextureBackground(tempTextureImageWindow);
}


void ImageWindow::CreateWindow()
{
#ifdef AMOVERBOSE
	std::cout<<"ImageWindow::CreateWindow()"<<std::endl;
#endif
	glutInitWindowSize(width, height);	
	window = glutCreateWindow(title.c_str());

}



	