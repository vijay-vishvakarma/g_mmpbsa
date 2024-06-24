#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# # Set up the module library
# add_library(pbcutil INTERFACE)

# set(PBCUTIL_SOURCES
#     ${CMAKE_CURRENT_SOURCE_DIR}/boxutilities.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/com.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/mshift.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/pbc.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/pbc_simd.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/pbcenums.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/pbcmethods.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/rmpbc.cpp
# )

# set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${PBCUTIL_SOURCES} PARENT_SCOPE)

# # Set GROMACS include directories
# find_package(GROMACS REQUIRED)
# target_include_directories(pbcutil INTERFACE
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     ${GROMACS_INCLUDE_DIRS}
# )

# # Link against GROMACS libraries
# target_link_libraries(pbcutil INTERFACE
#     legacy_api
#     ${GROMACS_LIBRARIES}
# )

# if (BUILD_TESTING)
#     add_subdirectory(tests)
# endif()
cmake_minimum_required(VERSION 3.29)

project(g_mmpbsa)

set(g_mmpbsa_MAJOR_VERSION 1)
set(g_mmpbsa_MINOR_VERSION 7)
set(g_mmpbsa_PATCH_VERSION 1)
set(g_mmpbsa_VERSION ${g_mmpbsa_MAJOR_VERSION}.${g_mmpbsa_MINOR_VERSION}.${g_mmpbsa_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
    set(DEF_INSTALL_CMAKE_DIR CMake)
else()
    set(DEF_INSTALL_CMAKE_DIR lib/CMake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
    set(var INSTALL_${p}_DIR)
    if(NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
endforeach()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set up include-directories
set(BUILD_INCLUDE ${CMAKE_BINARY_DIR}/include)
include_directories(
    "${PROJECT_SOURCE_DIR}"   # to find src/*.h
    "${PROJECT_BINARY_DIR}"   # to find src/config.h
    "${BUILD_INCLUDE}"        # to find header fix-ups
    "${CMAKE_CURRENT_SOURCE_DIR}/include" # include directory for pbcutil
)

# Specify the path to the GROMACS installation
set(GROMACS_DIR /opt/gromacs/installation)

# Find GROMACS package
find_package(GROMACS REQUIRED)

# Add sub-directories
add_subdirectory(src)
add_subdirectory(pbcutil)  # Add pbcutil library

# Link pbcutil and GROMACS libraries to the main executable
target_link_libraries(g_mmpbsa PRIVATE pbcutil ${GROMACS_LIBRARIES})

# ============== Eclipse CDT Error Handling======
#IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCC)
#IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
# ===============================
