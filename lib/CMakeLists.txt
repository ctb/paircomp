MESSAGE(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# these source files contain the templated version of paircomp
SET(SRCS
  MutableComparison.cc 
  ImmutableComparison.cc
  NwayComparison.cc
  algorithms.cc
  algorithms2.cc
)

ADD_LIBRARY(paircomplib STATIC ${SRCS})

#IF(APPLE)
#  SET_SOURCE_FILES_PROPERTIES(${SRCS}
#    PROPERTIES
#    COMPILE_FLAGS "-isysroot /Developer/SDKS/MacOSX10.4u.sdk -arch ppc -arch i386"
#    LINK_FLAGS "-Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch i386"
#  )
#ENDIF(APPLE)

ADD_SUBDIRECTORY( c++-tests )

