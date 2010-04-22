#  Copyright Haruhiko Matsuo and Synsge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps-looper)

######################################################################
# options
######################################################################

option(LOOPER_BUILD_TESTS "Build looper tests" ON)
option(LOOPER_BUILD_EXTRAS "Build in extras subdirectory" ON)

option(LOOPER_ENABLE_DIAG "Enable exact diagonalization method" ON)
option(LOOPER_ENABLE_ISING "Enable Swendsen-Wang method for classical Ising model" ON)
option(LOOPER_ENABLE_PI "Enable path integral method" ON)
option(LOOPER_ENABLE_SSE "Enable SSE method" ON)
option(LOOPER_ENABLE_SSE_QWL "Enable SSE quantum Wang-Landau method" OFF)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/ALPS" CACHE PATH "ALPS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

######################################################################
# Version information
######################################################################

set(LOOPER_VERSION_MAJOR 3)
set(LOOPER_VERSION_MINOR 2b14)
set(LOOPER_VERSION_BUILD "")

if(NOT LOOPER_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
  find_package(Subversion) 
  if(Subversion_FOUND) 
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} LOOPER) 
    set(LOOPER_VERSION_BUILD "r${LOOPER_WC_REVISION}")
  endif(Subversion_FOUND) 
endif(NOT LOOPER_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

if(LOOPER_VERSION_BUILD)
  set(LOOPER_VERSION "${LOOPER_VERSION_MAJOR}.${LOOPER_VERSION_MINOR}-${LOOPER_VERSION_BUILD}")
else(LOOPER_VERSION_BUILD)
  set(LOOPER_VERSION "${LOOPER_VERSION_MAJOR}.${LOOPER_VERSION_MINOR}")
endif(LOOPER_VERSION_BUILD)
message(STATUS "Looper version: ${LOOPER_VERSION}")

configure_file(looper/version.h.in ${PROJECT_BINARY_DIR}/looper/version.h)

######################################################################
# find ALPS Library
######################################################################

find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "ALPS version: ${ALPS_VERSION}")
include(${ALPS_USE_FILE})
if(NOT LAPACK_FOUND)
  message(FATAL "loop will not be built since lapack library is not found")
endif(NOT LAPACK_FOUND)

######################################################################
# applications
######################################################################

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

set(LOOP_CORE_SRC loop_custom.C loop_model.C)
if(LOOPER_ENABLE_DIAG)
  set(LOOP_CORE_SRC ${LOOP_CORE_SRC} diag.C)
endif(LOOPER_ENABLE_DIAG)
if(LOOPER_ENABLE_ISING)
  set(LOOP_CORE_SRC ${LOOP_CORE_SRC} ising.C)
endif(LOOPER_ENABLE_ISING)
if(LOOPER_ENABLE_PI)
  set(LOOP_CORE_SRC ${LOOP_CORE_SRC} path_integral.C path_integral_omp.C loop_model_omp.C)
endif(LOOPER_ENABLE_PI)
if(LOOPER_ENABLE_SSE)
  set(LOOP_CORE_SRC ${LOOP_CORE_SRC} sse.C)
endif(LOOPER_ENABLE_SSE)
if(LOOPER_ENABLE_SSE_QWL)
  set(LOOP_CORE_SRC ${LOOP_CORE_SRC} sse_qwl.C)
endif(LOOPER_ENABLE_SSE_QWL)

add_library(loop_core ${LOOP_CORE_SRC})
target_link_libraries(loop_core ${ALPS_LIBRARIES})

add_executable(loop loop.C)
target_link_libraries(loop loop_core)

######################################################################
# test
######################################################################

enable_testing()

if(LOOPER_BUILD_TESTS)
  set(TESTS loop_ns)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.C)
    target_link_libraries(${name} loop_core)
    add_test(${name}
      ${CMAKE_COMMAND}
        -Dcmd=${name}
        -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
        -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
        -Ddllexedir=${PROJECT_BINARY_DIR}/bin
        -Dinput=${name}.ip
        -Doutput=${name}.op
        -Dtimeout=600
        -P ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake
      )
  endforeach(name)
endif(LOOPER_BUILD_TESTS)

if(LOOPER_BUILD_TESTS)
  add_subdirectory(test)
endif(LOOPER_BUILD_TESTS)

######################################################################
# extras
######################################################################

if(LOOPER_BUILD_EXTRAS)
  add_subdirectory(extras)
endif(LOOPER_BUILD_EXTRAS)

######################################################################
# standalone
######################################################################

# Standalone version is built completeley independetly.
# add_subdirectory(standalone)
