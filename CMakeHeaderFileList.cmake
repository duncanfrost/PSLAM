# Set HEADER_subdir variable to all the files we want to parse during
# the build process. 
# Don't forget to update HEADER_ALL variable if you add/remove a 
# HEADER_subdir variable
#
# If you add/remove a directory, modify here

set(HEADER_OBJ
	AmoDefines.h
	AmoDefinesParameters.h
	Runnable.h
	CvWrapper.h
        AmoTimer.h
)


set(HEADER_MAPENGINES
	MapEngines/HomogFinder.h
	MapEngines/GeoFunctions.h
	MapEngines/MapOptimiser.h
	MapEngines/BundleAdjuster.h
	MapEngines/MapOptimiserEssential.h
	MapEngines/GraphFunctions.h
	MapEngines/PoseGraphOpt.h
)
set(HEADER_TRACKENGINES
	TrackEngines/ImageProcess.h
	TrackEngines/GlobalTransfoEstim.h
	TrackEngines/MapTracker.h
	TrackEngines/RobustMatching.h
)
set(HEADER_PRIMITIVES
	Primitives/Camera.h
	Primitives/HomogeneousMatrix.h
	Primitives/KeyFrame.h
	Primitives/obMap.h
	Primitives/MapPoint.h
)
set(HEADER_VISU
	Visualisation/VisualisationModule.h
	Visualisation/MotherWindow.h
	Visualisation/EmptyWindow.h
	Visualisation/ImageWindow.h
	Visualisation/MapWindow.h
	Visualisation/glTools.h
	Visualisation/glTexture.h
	Visualisation/objloader.h
	Visualisation/PlottingWindow.h
	Visualisation/VisuButton.h
	Visualisation/SimuWindow.h
)
set(HEADER_SRC
	ImageSource/VideoSource.h  
	ImageSource/VideoSourceLiveCV.h  
	ImageSource/VideoSourceSeq.h  
)

set(HEADER_VISO
	VisuOdo/filter.h   
	VisuOdo/matrix.h  
	VisuOdo/timer.h  
	VisuOdo/triangle.h
  	VisuOdo/viso.h
	VisuOdo/matcher.h
	VisuOdo/motionEstim.h
	VisuOdo/matcherWrapper.h
)


SET (HEADER_ALL 
  ${HEADER_OBJ}
  ${HEADER_MAPENGINES}
  ${HEADER_TRACKENGINES}
  ${HEADER_PRIMITIVES}
  ${HEADER_VISU}
  ${HEADER_SRC}
  ${HEADER_VISO}
  )
