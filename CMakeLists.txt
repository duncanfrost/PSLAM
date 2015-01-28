#define CMAKE_CXX_FLAGS_RELEASE as -O3 -DNDEBUG -fopenmp  for openmp use
#cmake_minimum_required(VERSION 2.5)
#put #define EIGEN_DONT_ALIGN in Eigen/src/Core/util/Macros.h

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


PROJECT(OBJECTR3D)
INCLUDE (FindPkgConfig)
SET (CMAKE_VERBOSE_MAKEFILE true)


SET(COMPUTEROFFICE 1)

#-----------------------------------------------------------------------------
# OBJECTR3D version number.  An even minor number corresponds to releases.
SET(OBJECTR3D_VERSION_MAJOR 1)
SET(OBJECTR3D_VERSION_MINOR 0)
SET(OBJECTR3D_VERSION_PATCH 0)
SET(OBJECTR3D_VERSION "${OBJECTR3D_VERSION_MAJOR}.${OBJECTR3D_VERSION_MINOR}")
SET(OBJECTR3D_VERSION_FULL "${OBJECTR3D_VERSION}.${OBJECTR3D_VERSION_PATCH}")      

SET(OBJECTR3D_INCLUDE_DIR "${OBJECTR3D_BINARY_DIR}/include/tracking")
SET(OBJECTR3D_DOC_DIR "${OBJECTR3D_BINARY_DIR}/doc")

#-----------------------------------------------------------------------------

# The location in which to install Example libraries.
SET(LIBRARY_OUTPUT_PATH ${OBJECTR3D_BINARY_DIR}/lib)

# The location in which to install OBJECTR3D binaries.
SET(BINARY_OUTPUT_PATH ${OBJECTR3D_BINARY_DIR}/bin)

# the include directory we depend on
SET(OBJECTR3D_INTERN_INCLUDE_DIR ${OBJECTR3D_BINARY_DIR}/include)
# the third party include directories we depend on
SET(OBJECTR3D_EXTERN_INCLUDE_DIR "")
# the internal visp libs we want to createn
SET(OBJECTR3D_INTERN_LIBS tracking)

# the third party libraries we depend on
SET(OBJECTR3D_EXTERN_LIBS "")
# the third party libraries we depend on
SET(OBJECTR3D_DEFS "")

# Get the OS
SET(OS ${CMAKE_SYSTEM_NAME})

## where are user-specific cmake modules
## make_dist.sh relies on $ENV{CMAKE_MODULES_PATH} (TODO, FIXME)
SET(CMAKE_MODULE_PATH ${OBJECTR3D_SOURCE_DIR}/CMakeModules)

INCLUDE(CheckLibraryExists)

# Create include directory which will contain all the headers
MAKE_DIRECTORY(${OBJECTR3D_INCLUDE_DIR})



#---------------------------------------------------------------------
# Platform dependent 
#----------------------------------------------------------------------
IF(WIN32 AND NOT CYGWIN)
  SET(OBJECTR3D_DEFS "${OBJECTR3D_DEFS} -DWIN32")
  SET(OBJECTR3D_DEFS "${OBJECTR3D_DEFS} -DGX") # To support exceptions
  CHECK_LIBRARY_EXISTS("winmm"  getch "" HAVE_LIBWINMM) # for timeGetTime()
  IF(HAVE_LIBWINMM)
    SET(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} "winmm")
  ENDIF(HAVE_LIBWINMM)
  IF(MSVC) #check if MSVC is the compiler and set ART_USE_MSVC accordingly
    SET(ART_USE_MSVC   TRUE)  # for header mbtConfig.h
  ENDIF(MSVC)
 IF(MSVC) #check if MSVC is the compiler and set ART_USE_MSVC accordingly
    SET(MBT_USE_MSVC   TRUE)  # for header mbtConfig.h
  ENDIF(MSVC)
ENDIF(WIN32 AND NOT CYGWIN)
IF(UNIX)
  SET(OBJECTR3D_DEFS "${OBJECTR3D_DEFS} -DUNIX")
  IF(APPLE)
    # We will not use rpath support.  Tell CMake to skip it.
    SET(CMAKE_SKIP_RPATH 1 CACHE INTERNAL "Whether to build with rpath." FORCE)
    SET(OBJECTR3D_DEFS "${OBJECTR3D_DEFS} -DAPPLE")
  ENDIF(APPLE)
ENDIF(UNIX)

#----------------------------------------------------------------------
# Third party libraries
#----------------------------------------------------------------------
# use sse3 instruction set
#SET(CMAKE_CXX_FLAGS "-msse3")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fopenmp -msse3")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -DNDEBUG")

#add_definitions(USE_CVD)
#set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} -lcvd)

#to remove warning deprecated from interpolation:
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-write-strings")

#include eigen
include_directories("/home/amaury/Program/Libs/eigen-eigen-ffa86ffb5570")
#include cuda just for fevtor uchar3 definition
include_directories("/usr/local/cuda-5.5/targets/x86_64-linux/include/")



#----------------------------------------------------------------------
# GL
#----------------------------------------------------------------------
find_package(GLUT)
find_package(OpenGL)
set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} ${GLUT_LIBRARY} ${OPENGL_LIBRARY})
set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} -lGLEW)

#----------------------------------------------------------------------
# OpenCV - For image acquisition
#----------------------------------------------------------------------

pkg_search_module (OPENCV opencv REQUIRED)
include_directories(${OPENCV_INCLUDE_DIRS})
set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} ${OPENCV_LDFLAGS})



message( "${OPENCV_LDFLAGS}")

#set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} -lopencv_core -lopencv_ml -lopencv_highgui -lopencv_imgproc -lopencv_features2d -lopencv_gpu -lopencv_flann -lopencv_nonfree -lopencv_ocl -lopencv_calib3d -lopencv_video)



#if opencv searched for cuda nppi
LINK_DIRECTORIES( "/usr/local/cuda/targets/x86_64-linux/lib/" )
include_directories("/usr/local/cuda/targets/x86_64-linux/include")




#----------------------------------------------------------------------
# Third party libraries
#----------------------------------------------------------------------


#------------------------------------------------------------
# OpenMP
#------------------------------------------------------------
IF(CUDA_FOUND)
OPTION(USE_OMP "Use OpenMP" off)
ELSE(CUDA_FOUND)
OPTION(USE_OMP "Use OpenMP" ON)
ENDIF(CUDA_FOUND)

#set(OBJECTR3D_EXTERN_LIBS ${OBJECTR3D_EXTERN_LIBS} -fopenmp)

OPTION(USE_ICC "Use Intel Compiler" ON)

IF(USE_OMP)
    IF(UNIX)
        IF(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dnoomp -openmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dnoomp -openmp")
        ELSE(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        ENDIF(USE_ICC)
    ELSE(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
    ENDIF(UNIX)
    ADD_DEFINITIONS(-D USE_OMP)
ENDIF(USE_OMP)

#----------------------------------------------------------------------
# Concatenate intern and third party include dirs
#----------------------------------------------------------------------
INCLUDE_DIRECTORIES(${OBJECTR3D_INTERN_INCLUDE_DIR} ${OBJECTR3D_EXTERN_INCLUDE_DIR})

#----------------------------------------------------------------------
# Add library path
#----------------------------------------------------------------------
LINK_DIRECTORIES(${OBJECTR3D_SOURCE_DIR}/lib)

#----------------------------------------------------------------------
# Add definitions
#----------------------------------------------------------------------
ADD_DEFINITIONS(${OBJECTR3D_DEFS})

IF ($ENV{TRAVAIL_LOCAL})
   ADD_DEFINITIONS(-D USE_TRAVAIL_LOCAL)
ENDIF($ENV{TRAVAIL_LOCAL})

IF ($ENV{TRAVAIL_LOCAL_DIONE})
   ADD_DEFINITIONS(-D USE_TRAVAIL_LOCAL_DIONE)
ENDIF($ENV{TRAVAIL_LOCAL_DIONE})

#----------------------------------------------------------------------
# Propagation in sub dirs
#----------------------------------------------------------------------
SUBDIRS(src)
SUBDIRS(BinSamples)
SUBDIRS(BinSlam)
SUBDIRS(BinTest)



