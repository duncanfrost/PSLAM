# Set SRC_subdir variable to all the files we want to parse during
# the build process. 
# Don't forget to update SRC_ALL variable if you add/remove a SRC_subdir 
# variable
#
# If you add/remove a directory, modify here


SET(SRC_OBJ
	AmoDefines.cpp
	CvWrapper.cpp
)

set(SRC_MAPENGINES
	MapEngines/HomogFinder.cpp
	MapEngines/GeoFunctions.cpp
	MapEngines/MapOptimiser.cpp
	MapEngines/BundleAdjuster.cpp
	MapEngines/MapOptimiserEssential.cpp
	MapEngines/PoseGraphOpt.cpp
)
set(SRC_TRACKENGINES
	TrackEngines/ImageProcess.cpp
	TrackEngines/GlobalTransfoEstim.cpp
	TrackEngines/MapTracker.cpp
	TrackEngines/RobustMatching.cpp
)
set(SRC_PRIMITIVES
	Primitives/Camera.cpp
	Primitives/HomogeneousMatrix.cpp
	Primitives/KeyFrame.cpp
	Primitives/obMap.cpp
	Primitives/MapPoint.cpp
)
set(SRC_VISU
	Visualisation/VisualisationModule.cpp
	Visualisation/MotherWindow.cpp
	Visualisation/EmptyWindow.cpp
	Visualisation/ImageWindow.cpp
	Visualisation/MapWindow.cpp
	Visualisation/SimuWindow.cpp
	Visualisation/glTools.cpp
	Visualisation/glTexture.cpp
	Visualisation/objloader.cpp
	Visualisation/PlottingWindow.cpp
 	Visualisation/VisuButton.cpp
)
SET(SRC_SRC
	ImageSource/VideoSource.cpp
	ImageSource/VideoSourceLiveCV.cpp  
	ImageSource/VideoSourceSeq.cpp  
)
set(SRC_VISO
	VisuOdo/filter.cpp  
	VisuOdo/matcher.cpp 
	VisuOdo/matrix.cpp 
	VisuOdo/triangle.cpp
  	VisuOdo/viso.cpp
  	VisuOdo/motionEstim.cpp
	VisuOdo/matcherWrapper.cpp
)

SET (SRC_ALL
  ${SRC_OBJ}
  ${SRC_MAPENGINES}
  ${SRC_TRACKENGINES}
  ${SRC_PRIMITIVES}
  ${SRC_VISU}
  ${SRC_SRC}
  ${SRC_VISO}
  )
