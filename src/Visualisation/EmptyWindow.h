#pragma once

/*#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/rgba.h>
*/
#include <objectr3d/MotherWindow.h>

//TextureSet imgTexture;
	
class EmptyWindow:public MotherWindow
{
 public:
	EmptyWindow(std::string _title,int _width,int _height);
	void CreateWindow();
	void prepare_draw(){};
	
 protected: 
	void (*DrawFunction)();
};
