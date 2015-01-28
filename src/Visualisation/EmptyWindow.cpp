
#include "EmptyWindow.h"

EmptyWindow::EmptyWindow(std::string _title,int _width,int _height)
{
	title=_title;
	width=_width;
	height=_height;
}

void EmptyWindow::CreateWindow()
{
#ifdef AMOVERBOSE
	std::cout<<"EmptyWindow::CreateWindow()"<<std::endl;
#endif
	glutInitWindowSize(width, height);	
	window = glutCreateWindow(title.c_str());

}



	