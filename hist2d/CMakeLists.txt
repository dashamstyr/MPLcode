project(fasthist)
cmake_minimum_required(VERSION 2.8)
message("system is ${CMAKE_SYSTEM_NAME}")
set(homedir $ENV{HOME})
set(CMAKE_MODULE_PATH "${homedir}/cmake")
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/testlib)
set(CMAKE_INSTALL_PREFIX ${homedir}/lib/python)
set_source_files_properties(
  ${fasthist_SOURCE_DIR}/fasthist.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )


include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH} ${CMAKE_SOURCE_DIR})
include(UseCython)
message("python is ${PYTHON_EXECUTABLE}")


set(CYTHON_CXX_EXTENSION "cxx")

find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)
find_package(PyInclude REQUIRED)
include_directories(${TRUE_PYTHON_INC} ${NUMPY_INCLUDE_PATH} ${fasthist_SOURCE_DIR}})


cython_add_module(fasthist fasthist.pyx)

set_target_properties(fasthist PROPERTIES PREFIX "")
install(TARGETS fasthist LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})

add_custom_target(test COMMAND ${PYTHON_EXECUTABLE} ${fasthist_SOURCE_DIR}/test_hist.py 
                    fasthist ${LIBRARY_OUTPUT_PATH}
                    DEPENDS fasthist COMMENT "execute test_hist.py")
