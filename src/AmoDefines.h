#pragma once

//#define AMOVERBOSE

/*#ifdef USE_FREEGLUT
	#define HGLUT <GL/freeglut.h>
#else*/
	#define HGLUT <GL/glut.h>
//#endif

//#define EIGEN_DONT_ALIGN

#define WORK_DESKTOP
#ifdef WORK_DESKTOP
  #define PEYE_CALIB_FILE "camera.cfg"
  #define KITTI_CALIB_FILE "camera2.cfg"
  #define WIMU_WORLD_DIR "/home/amaury/Program/Projects/files/Simu/noTransparent/"
#else
  #define PEYE_CALIB_FILE "/home/momo/Progam/Data/camera.cfg"
  #define WIMU_WORLD_DIR "/home/momo/Progam/Data/noTransparent/"
#endif


#include "AmoDefinesParameters.h"

#define OUTLIERCOUNTMAX 5


#define SAVE_POINT_COLOR

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "AmoTimer.h"

void InitProcAndGPU();

#include <iostream>
#include <iomanip>

#define coutTracking std::cout <<std::right << std::setw(90)
#define endlTracking std::endl
#define coutColTracking std::cout << "\033[1;30m" <<std::right << std::setw(90)
#define endlColTracking "\033[0m"<<std::endl
#define coutErrTracking std::cout << "\033[1;31m" <<std::right << std::setw(90)
#define endlErrTracking "\033[0m"<<std::endl


#define coutMapping std::cout
#define endlMapping std::endl
#define coutColMapping std::cout << "\033[1;32m"
#define endlColMapping "\033[0m"<<std::endl
#define coutErrMapping std::cout << "\033[1;31m"
#define endlErrMapping "\033[0m"<<std::endl

#define coutGreen std::cout << "\033[1;32m"
#define endlGreen "\033[0m"<<std::endl
#define coutRed std::cout << "\033[1;31m"
#define endlRed "\033[0m"<<std::endl


