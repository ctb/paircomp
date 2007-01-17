INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

ADD_EXECUTABLE(paircomp paircomp.cc)
TARGET_LINK_LIBRARIES(paircomp paircomplib)

ADD_EXECUTABLE(seqcomp seqcomp.cc)
TARGET_LINK_LIBRARIES(seqcomp paircomplib)

ADD_EXECUTABLE(find_patch find_patch.cc)
TARGET_LINK_LIBRARIES(find_patch paircomplib)

ADD_EXECUTABLE(window-extract window-extract.cc)
TARGET_LINK_LIBRARIES(window-extract paircomplib)

INSTALL(TARGETS paircomp seqcomp find_patch
        RUNTIME DESTINATION bin
)
