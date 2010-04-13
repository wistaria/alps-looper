cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(alps-looper)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/ALPS" CACHE PATH "ALPS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/CMake)

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
include_directories(${PROJECT_SOURCE_DIR})

add_library(loop_common loop_custom.C loop_model.C diag.C ising.C path_integral.C sse.C)
target_link_libraries(loop_common ${ALPS_LIBRARIES})

add_executable(loop loop.C)
target_link_libraries(loop loop_common)
install(TARGETS loop RUNTIME DESTINATION bin)

######################################################################
# tests
######################################################################

enable_testing()
add_custom_target(check ${CMAKE_CTEST_COMMAND})
add_subdirectory(test EXCLUDE_FROM_ALL)

set(TESTS loop_ns)
foreach(name ${TESTS})
  add_executable(${name} EXCLUDE_FROM_ALL ${name}.C)
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
  add_dependencies(check ${name})
endforeach(name)

######################################################################
# subdirectories
######################################################################

