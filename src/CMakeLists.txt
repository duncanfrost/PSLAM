#----------------------------------------------------------------------
# Include project files 
#----------------------------------------------------------------------
# include all the Example project .cpp files
INCLUDE(${OBJECTR3D_SOURCE_DIR}/CMakeSourceFileList.cmake)
# include all the Example project .h files
INCLUDE(${OBJECTR3D_SOURCE_DIR}/CMakeHeaderFileList.cmake)

#----------------------------------------------------------------------
# Create rule to copy all the headers from src to include/tracking
#----------------------------------------------------------------------
# For each header, we create a rule
SET(HEADER_IN_INCLUDE_DIR "")
SET(ADDED_SRC_DEPEND_PROPERTIES "")
FOREACH(header ${HEADER_ALL})
  GET_FILENAME_COMPONENT(headerName ${header} NAME) 
  ADD_CUSTOM_COMMAND(
    OUTPUT ${OBJECTR3D_BINARY_DIR}/include/objectr3d/${headerName}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${header} ${OBJECTR3D_BINARY_DIR}/include/objectr3d/${headerName}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${header}
    )
  SET(HEADER_IN_INCLUDE_DIR ${HEADER_IN_INCLUDE_DIR} ${OBJECTR3D_BINARY_DIR}/include/objectr3d/${headerName}
    )
ENDFOREACH(header)





#----------------------------------------------------------------------
# Create a top level rule to copy all the headers from src to include/tracking
#----------------------------------------------------------------------
ADD_CUSTOM_TARGET(header ALL
  DEPENDS ${HEADER_IN_INCLUDE_DIR}
  )


#----------------------------------------------------------------------
# build rule for the library
#----------------------------------------------------------------------
ADD_LIBRARY (${OBJECTR3D_INTERN_LIBS}
  ${SRC_ALL} 
  ${HEADER_ALL}
  ) 


# create the headers in include/tracking before compiling the lib
ADD_DEPENDENCIES(${OBJECTR3D_INTERN_LIBS} header)

# Append the library version information to the library target
# properties. 
SET_TARGET_PROPERTIES(${OBJECTR3D_INTERN_LIBS}
  PROPERTIES
  VERSION ${OBJECTR3D_VERSION_MAJOR}.${OBJECTR3D_VERSION_MINOR}
  SOVERSION ${OBJECTR3D_VERSION_PATCH}
  )

#----------------------------------------------------------------------
# customize install target 
#----------------------------------------------------------------------
# install rule for the library
INSTALL_TARGETS(/lib ${OBJECTR3D_INTERN_LIBS})
# install rule for all the headers
INSTALL_FILES(/include/tracking FILES ${HEADER_ALL})
