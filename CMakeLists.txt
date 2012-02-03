cmake_minimum_required(VERSION 2.6)

project(DIPBAC)
# Name of the project this creates DIPBAC_SOURCE_DIR as being the location
# of this file

# Built-in modules
include(CheckIncludeFiles)
include(CheckLibraryExists)

# Set up the directory strucure
set(LOGDIRECTORY "${DIPBAC_SOURCE_DIR}/log")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${DIPBAC_SOURCE_DIR}/bin")

# Add extra directories for finding header files
include_directories("${DIPBAC_SOURCE_DIR}/include")
include_directories("${DIPBAC_SOURCE_DIR}/CPP_Base/include")
include_directories("/usr/local/include/eigen3")

# Turn on debugging warnings and profiling
if(CMAKE_BUILD_TYPE MATCHES "Develop")
    set(CMAKE_CXX_FLAGS_DEVELOP "${CMAKE_CXX_FLAGS}")

    # Build a debugging version
    add_executable(DiPBaCpp_Debug DiPBaC.cpp)
    set_target_properties(DiPBaCpp_Debug PROPERTIES COMPILE_FLAGS "-g -Wall" LINK_FLAGS "-g -Wall")	
    
    # Build a profiling version
    add_executable(DiPBaCpp_Profile DiPBaC.cpp)
    set_target_properties(DiPBaCpp_Profile PROPERTIES COMPILE_FLAGS "-pg -g -Wall"  LINK_FLAGS "-pg -g -Wall")	
    
endif(CMAKE_BUILD_TYPE MATCHES "Develop")

# Build as optimised code
add_executable(DiPBaCpp DiPBaC.cpp)
set_target_properties(DiPBaCpp PROPERTIES COMPILE_FLAGS "-O3"  LINK_FLAGS "-O3")	


