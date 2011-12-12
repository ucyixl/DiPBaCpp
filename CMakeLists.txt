cmake_minimum_required(VERSION 2.6)

project(PROFREGR)
# Name of the project this creates PROFREGR_SOURCE_DIR as being the location
# of this file

# Turn on debugging warnings and profiling
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -Wall")
# Build as optimised code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Built-in modules
include(CheckIncludeFiles)
include(CheckLibraryExists)

# Set up the directory strucure
set(LOGDIRECTORY "${PROFREGR_SOURCE_DIR}/log")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROFREGR_SOURCE_DIR}/modules/")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROFREGR_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROFREGR_SOURCE_DIR}/bin")

# Make sure we build shared libraries
set(BUILD_SHARED_LIBS ON)

# Add extra directories for finding header files
include_directories("${PROFREGR_SOURCE_DIR}/include")

# Add extra directories for finding libraries
link_directories("${PROFREGR_SOURCE_DIR}/lib")

# Testing
#enable_testing()
# This command means "make check" will build the tests 
# "make test" will be used to run the test
#add_custom_target(check COMMAND ${CMAKE_TEST_COMMAND})

add_executable(profileRegression profileRegression.cpp)
#set_target_properties(profileRegression PROPERTIES COMPILE_FLAGS "-fopenmp")
target_link_libraries (profileRegression armadillo)
target_link_libraries (profileRegression goto2)
target_link_libraries (profileRegression gfortran)
