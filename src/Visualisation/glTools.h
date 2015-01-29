/*This file defines several functions used to simultaneously use
 * Toon and Opengl, a structure that regroups all the information
 * of the camera (intrinsic and extrinsic params) in Toon and GL
 * style. 
 * It also defines a class Camera used by the viewers and a structure
 * to save textures to project them on surfaces.*/

#pragma once

#include "../AmoDefines.h"
#include "../Primitives/Camera.h"
#include "../Primitives/HomogeneousMatrix.h"
#include HGLUT

#include <iostream>
#include <fstream>
#include <memory.h>
#include <cmath>

#include <opencv/cv.h>
#include <opencv/highgui.h>

void multiplyGlMatrix(Matrix4f _mat);

//set projection matrix and modelview of opengl from cam, to use primitives as glVertex3f
void set3DGLProjection(Camera *_cam,HomogeneousMatrix _modelView);
void unset3DGLProjection();


//set 2D projection matrix to use things as glVertex2f (top left=0,0)
void set2DGLProjection();
void unset2DGLProjection();

//screenshot functions
cv::Mat Capture_Image();
cv::Mat Capture_ZBuffer();
cv::Mat Capture_realZ();
cv::Mat Capture_displayable_ZBuffer();

float zbuffer_to_zreal(float _z);
float zreal_to_zbuffer(float _z);

void floatImageSave(const char* szFilename, cv::Mat &img);
bool floatImageRead(const char* szFilename, cv::Mat &img);

void getDisplayableImage(cv::Mat &img,cv::Mat &img_disp,float ignored_value=-1e10);

void drawText(GLint x, GLint y, const char* s, GLfloat r, GLfloat g, GLfloat b);
