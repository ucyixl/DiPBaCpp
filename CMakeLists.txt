# (C) Copyright David Hastie and Silvia Liverani, 2012.

# DiPBaC++ is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.

# DiPBaC++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with DiPBaC++ in the documentation directory. If not, see
# <http://www.gnu.org/licenses/>.

# The external linear algebra library Eigen, parts of which are included  in the
# lib directory is released under the LGPL3+ licence. See comments in file headers
# for details.

# The Boost C++ header library, parts of which are included in the  lib directory
# is released under the Boost Software Licence, Version 1.0, a copy  of which is
# included in the documentation directory.

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

# Extract the boost files
execute_process(
COMMAND ${CMAKE_COMMAND} -E tar xzf "${DIPBAC_SOURCE_DIR}/lib/include/boost_1.48.0.tar.gz"
WORKING_DIRECTORY "${DIPBAC_SOURCE_DIR}/lib/include"
)

# Extract the Eigen files
execute_process(
COMMAND ${CMAKE_COMMAND} -E tar xzf "${DIPBAC_SOURCE_DIR}/lib/include/Eigen_3.0.5.tar.gz"
WORKING_DIRECTORY "${DIPBAC_SOURCE_DIR}/lib/include"
)


# Add extra directories for finding header files
# The before key word means that the libraries include in this package are 
# used instead of any already installed on the system 
include_directories(BEFORE "${DIPBAC_SOURCE_DIR}/lib/include")
include_directories("${DIPBAC_SOURCE_DIR}/include")

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


