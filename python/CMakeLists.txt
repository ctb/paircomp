FIND_PACKAGE(PythonInterp)

IF(PYTHONINTERP_FOUND) 
  # I needed to split into the C compiling and python compiling
  # so I could add the -L option to specify where CMake had built the
  # library (-L option)
  ADD_CUSTOM_TARGET(paircomp_ext ALL
       COMMAND ${PYTHON_EXECUTABLE} 
               ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
               build_ext
               -t ${CMAKE_CURRENT_BINARY_DIR}
               -b ${CMAKE_CURRENT_BINARY_DIR}
               -L ${CMAKE_CURRENT_BINARY_DIR}/../lib
               -f
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
       DEPENDS paircomplib
       VERBATIM
  )
  ADD_CUSTOM_TARGET(paircomp_py ALL
       COMMAND ${PYTHON_EXECUTABLE} 
               ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
               build_py
               -d ${CMAKE_CURRENT_BINARY_DIR}
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
       VERBATIM
  )
ENDIF(PYTHONINTERP_FOUND)
