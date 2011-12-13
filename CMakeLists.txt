cmake_minimum_required(VERSION 2.6)

project(DIPBAC)
# Name of the project this creates DIPBAC_SOURCE_DIR as being the location
# of this file

# Turn on debugging warnings and profiling
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -Wall")
# Build as optimised code
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Built-in modules
include(CheckIncludeFiles)
include(CheckLibraryExists)

# Set up the directory strucure
set(LOGDIRECTORY "${DIPBAC_SOURCE_DIR}/log")
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${DIPBAC_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${DIPBAC_SOURCE_DIR}/bin")

# Make sure we build shared libraries
#set(BUILD_SHARED_LIBS ON)

# Add extra directories for finding header files
include_directories("${DIPBAC_SOURCE_DIR}/include")
include_directories("${DIPBAC_SOURCE_DIR}/CPP_Base/include")

# Add extra directories for finding libraries
#link_directories("${DIPBAC_SOURCE_DIR}/lib")

# Testing
#enable_testing()
# This command means "make check" will build the tests 
# "make test" will be used to run the test
#add_custom_target(check COMMAND ${CMAKE_TEST_COMMAND})

add_executable(DiPBaCpp DiPBaC.cpp)
#set_target_properties(profileRegression PROPERTIES COMPILE_FLAGS "-fopenmp")
target_link_libraries (DiPBaCpp armadillo)
target_link_libraries (DiPBaCpp goto2)
target_link_libraries (DiPBaCpp gfortran)

