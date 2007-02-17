INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

# This test requires hard coded filesystem paths, which doesn't currently work
# with the CMake out-of-source builds
# ADD_EXECUTABLE(test-empty test-empty.cc)
# TARGET_LINK_LIBRARIES(test-empty paircomplib)               
# ADD_TEST(test-empty test-empty)

ADD_EXECUTABLE(test-transit test-transit.cc)
TARGET_LINK_LIBRARIES(test-transit paircomplib)
ADD_TEST(test-transit test-transit)
