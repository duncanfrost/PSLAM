#pragma once

/*#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/rgba.h>
*/
#include <objectr3d/AmoDefines.h>
#include HGLUT
#include <objectr3d/glTools.h>
//TextureSet imgTexture;
	
class VisuButton
{
 public:
	//create button with title, onClick and pos in window
	VisuButton(std::string _caption,void (*_onClick)(),int _x,int _y, int _w,int _h);
	void checkForAtivity(int button, int state, int x, int y);
	void checkForAtivityMotion(int x,int y);

	void draw();
	
 protected: 
	void (*onClick)();
	
	std::string caption;
	int x; int y; 
	int w; int h;
	
	bool clicked;
	bool mouse_in;
	
};
