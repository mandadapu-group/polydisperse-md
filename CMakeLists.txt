CMAKE_MINIMUM_REQUIRED(VERSION 3.9 FATAL_ERROR)

# Name the plugin project
project(polymd LANGUAGES C CXX)

# Find the installed HOOMD
find_package(HOOMD REQUIRED)

message(STATUS "Found HOOMD ${HOOMD_VERSION}: ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}")

set(ENABLE_HIP ${ENABLE_GPU})

# Force installation to the HOOMD installation location
set(CMAKE_INSTALL_PREFIX ${HOOMD_INSTALL_PREFIX} CACHE PATH "Installation prefix" FORCE)

# add subdirectories
add_subdirectory(${PROJECT_NAME})

#add_subdirectory (CMake)

#include (HOOMDHIPSetup)
#include (HOOMDCUDASetup)
