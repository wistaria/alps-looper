#  Copyright Haruhiko Matsuo and Synsge Todo 2011-2018.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(alps-looper NONE)

######################################################################
# Version information
######################################################################

set(LOOPER_VERSION 4.1.0)
message(STATUS "Looper version: ${LOOPER_VERSION}")

configure_file(looper/version.h.in ${PROJECT_BINARY_DIR}/looper/version.h)

######################################################################
# options
######################################################################

option(LOOPER_BUILD_TESTS "Build looper tests" ON)

option(LOOPER_ENABLE_DIAG "Enable exact diagonalization method" ON)
option(LOOPER_ENABLE_ISING "Enable Swendsen-Wang method for classical Ising model" ON)
option(LOOPER_ENABLE_PI "Enable path integral method" ON)
option(LOOPER_ENABLE_SSE "Enable SSE method" ON)
option(LOOPER_ENABLE_SSE_QWL "Enable SSE quantum Wang-Landau method" OFF)

option(LOOPER_ENABLE_OPENMP "OpenMP is enabled by default. You can turn off OpenMP by -DLOOPER_ENABLE_OPENMP=OFF" ON)
if(LOOPER_ENABLE_OPENMP)
  add_definitions(-DLOOPER_ENABLE_OPENMP)
endif(LOOPER_ENABLE_OPENMP)

option(LOOPER_ENABLE_TIMER "You can turn on Timer by -DLOOPER_ENABLE_TIMER=ON" OFF)
if(LOOPER_ENABLE_TIMER)
  add_definitions(-DALPS_ENABLE_TIMER)
endif(LOOPER_ENABLE_TIMER)

if(LOOPER_TIMER)
  option(LOOPER_TIMER_TRACE "You can turn on debug trace output in timer by -DLOOPER_TIMER_TRACE=ON" OFF)
  if(LOOPER_TIMER_TRACE)
    add_definitions(-DALPS_ENABLE_TIMER_TRACE)
  endif(LOOPER_TIMER_TRACE)

  option(LOOPER_TIMER_DETAILED "You can turn on detailed timer report by -DLOOPER_TIMER_DETAILED=ON" OFF)
  if(LOOPER_TIMER_DETAILED)
    add_definitions(-DALPS_ENABLE_TIMER_DETAILED)
  endif(LOOPER_TIMER_DETAILED)

  option(LOOPER_TIMER_BARRIER "You can turn on barrier synchronization in timer by -DLOOPER_TIMER_BARRIER=ON" OFF)
  if(LOOPER_TIMER_BARRIER)
    add_definitions(-DALPS_ENABLE_TIMER_BARRIER)
  endif(LOOPER_TIMER_BARRIER)

  option(LOOPER_ENABLE_PA1 "You can turn on PA-info by -DLOOPER_ENABLE_PA=ON" OFF)
  option(LOOPER_ENABLE_PA2 "You can turn on old-PA-info by -DLOOPER_ENABLE_PA=ON" OFF)
  option(LOOPER_ENABLE_RANGEDSAMP "You can turn on ranged sampling by -DLOOPER_ENABLE_RANGEDSAMP=ON" OFF)
  if(LOOPER_ENABLE_PA1)
    add_definitions(-DALPS_ENABLE_PA1)
  endif(LOOPER_ENABLE_PA1)
  if(LOOPER_ENABLE_PA2)
    add_definitions(-DALPS_ENABLE_PA2)
  endif(LOOPER_ENABLE_PA2)
  if(LOOPER_ENABLE_RANGEDSAMP)
    add_definitions(-DALPS_ENABLE_RANGEDSAMP)
  endif(LOOPER_ENABLE_RANGEDSAMP)
  if(LOOPER_ENABLE_PA1 AND LOOPER_ENABLE_RANGEDSAMP) 
    message(FATAL_ERROR "Both of sampler and hardware counter can't be enabled.")
  endif(LOOPER_ENABLE_PA1 AND LOOPER_ENABLE_RANGEDSAMP) 
  if(LOOPER_ENABLE_PA2 AND LOOPER_ENABLE_RANGEDSAMP) 
    message(FATAL_ERROR "Both of sampler and hardware counter can't be enabled.")
  endif(LOOPER_ENABLE_PA2 AND LOOPER_ENABLE_RANGEDSAMP) 
endif(LOOPER_TIMER)

option(LOOPER_USE_DETERMINISTIC_UNIFY "Use deterministic unify for debug" OFF)
if(LOOPER_USE_DETERMINISTIC_UNIFY)
  add_definitions(-DLOOPER_USE_DETERMINISTIC_UNIFY)
endif(LOOPER_USE_DETERMINISTIC_UNIFY)

######################################################################
# find ALPS Library
######################################################################

find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} $ENV{ALPS_ROOT} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})
if(NOT LAPACK_FOUND)
  message(FATAL "loop will not be built since lapack library is not found")
endif(NOT LAPACK_FOUND)

# enable C++ compilers
enable_language(CXX)

find_package(MPI)

# Boost
find_package(Boost REQUIRED COMPONENTS mpi system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

######################################################################
# applications
######################################################################

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

set(LOOP_SRC loop.C loop_custom.C loop_model.C)

if(LOOPER_ENABLE_DIAG)
  set(LOOP_SRC diag.C ${LOOP_SRC})
endif(LOOPER_ENABLE_DIAG)
if(LOOPER_ENABLE_ISING)
  set(LOOP_SRC ising.C ${LOOP_SRC})
endif(LOOPER_ENABLE_ISING)
if(LOOPER_ENABLE_PI)
  set(LOOP_SRC path_integral.C ${LOOP_SRC})
  if(MPI_FOUND)
    set(LOOP_SRC path_integral_mpi.C ${LOOP_SRC})
  endif(MPI_FOUND)
endif(LOOPER_ENABLE_PI)
if(LOOPER_ENABLE_SSE)
  set(LOOP_SRC sse.C ${LOOP_SRC})
endif(LOOPER_ENABLE_SSE)
if(LOOPER_ENABLE_SSE_QWL)
  set(LOOP_SRC sse_qwl.C ${LOOP_SRC})
endif(LOOPER_ENABLE_SSE_QWL)

add_executable(loop ${LOOP_SRC})
target_link_libraries(loop alps Boost::mpi ${ALPS_EXTRA_LIBRARIES})
install(TARGETS loop RUNTIME DESTINATION bin)

######################################################################
# test
######################################################################

set(TESTS loop)
foreach(name ${TESTS})
  add_alps_test(${name})
endforeach(name)

if(LOOPER_BUILD_TESTS)
  add_subdirectory(test)
endif(LOOPER_BUILD_TESTS)

######################################################################
# extras
######################################################################

add_subdirectory(extras)

######################################################################
# standalone
######################################################################

# Standalone version is built completeley independetly.
# add_subdirectory(standalone)

######################################################################
# headers
######################################################################

install(DIRECTORY looper ${CMAKE_BINARY_DIR}/looper
        DESTINATION include
        PATTERN ".DS_Store" EXCLUDE
        PATTERN "*.in" EXCLUDE)

######################################################################
# kernels
######################################################################

install(FILES diag.C ising.C loop.C loop_custom.C loop_model.C path_integral.C path_integral_mpi.C qwl_evaluate.C sse.C sse_qwl.C DESTINATION share/looper/kernel)

######################################################################
# CMake configuration
######################################################################

configure_file(config/looper-config.cmake.in ${PROJECT_BINARY_DIR}/config/looper-config.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/config/looper-config.cmake DESTINATION share/looper)
