cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps-looper)

######################################################################
# options
######################################################################

option(LOOPER_INSTALL "Install loop application" OFF)
option(LOOPER_BUILD_TESTS "Build looper tests" ON)

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
    MESSAGE("Current revision is r${LOOPER_WC_REVISION} (${LOOPER_WC_LAST_CHANGED_DATE})") 
    set(LOOPER_VERSION_BUILD "r${LOOPER_WC_REVISION}")
  endif(Subversion_FOUND) 
endif(NOT LOOPER_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

if(LOOPER_VERSION_BUILD)
  set(LOOPER_VERSION "${LOOPER_VERSION_MAJOR}.${LOOPER_VERSION_MINOR}-${LOOPER_VERSION_BUILD}")
else(LOOPER_VERSION_BUILD)
  set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}")
endif(LOOPER_VERSION_BUILD)
set(LOOPER_VERSION_STRING "ALPS/looper version ${LOOPER_VERSION}")

configure_file(looper/version.h.in ${PROJECT_BINARY_DIR}/looper/version.h)

######################################################################
# find ALPS Library
######################################################################

find_package(ALPS REQUIRED)
include(${ALPS_USE_FILE})
if(NOT LAPACK_FOUND)
  message(FATAL "loop will not be built since lapack library is not found")
endif(NOT LAPACK_FOUND)

######################################################################
# RPATH setting
######################################################################

set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

######################################################################
# applications
######################################################################

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

add_library(loop_common loop_custom.C loop_model.C diag.C ising.C path_integral.C sse.C)
target_link_libraries(loop_common ${ALPS_LIBRARIES})

add_executable(loop loop.C)
target_link_libraries(loop loop_common)
if(LOOPER_INSTALL)
  install(TARGETS loop RUNTIME DESTINATION bin)
endif(LOOPER_INSTALL)

######################################################################
# tests
######################################################################

enable_testing()
if(LOOPER_BUILD_TESTS)
  add_subdirectory(test)

  set(TESTS loop_ns)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.C)
    target_link_libraries(${name} loop_common)
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
