//Copyright 2008 Isis Innovation Limited 
/*This file defines the function necessary to allocate
 textures on the GPU, takes care of the resizing with 
 power of two dimensions and uses CVD as input*/


#ifndef GLTEXTURE__H
#define GLTEXTURE__H

#include <objectr3d/AmoDefines.h>
#include HGLUT

#include <opencv/cv.h>
#include <opencv/highgui.h>


#include <iostream>
#include <fstream>
#include <memory.h>
#include "math.h"

//#include <cvd/utility.h>

struct TextureSet  
{
	TextureSet():isAllocated(false){};
	~TextureSet(){if(isAllocated)glDeleteTextures(1,&Texture);};
	GLuint Texture;
	float Text_wu;
	float Text_hu;
	bool isAllocated;
};

void LoadTexture(cv::Mat &Img,TextureSet &Textureset);
void DrawTextureBackground(TextureSet &Tex);

#endif
