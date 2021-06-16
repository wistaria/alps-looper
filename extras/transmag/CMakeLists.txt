add_definitions(-DLOOP_CONFIG_HEADER="extras/transmag/transmag_config.h")

set(LOOP_SRC ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/diag.C ${PROJECT_SOURCE_DIR}/path_integral.C)
if(MPI_FOUND)
  set(LOOP_SRC ${LOOP_SRC} ${PROJECT_SOURCE_DIR}/path_integral_mpi.C)
endif(MPI_FOUND)
if(CMAKE_ASM_COMPILER)
  enable_language(ASM)
  set_property(SOURCE ${PROJECT_SOURCE_DIR}/looper/atomic.s PROPERTY LANGUAGE ASM)
  set(LOOP_SRC ${LOOP_SRC} ${PROJECT_SOURCE_DIR}/looper/atomic.s)
endif(CMAKE_ASM_COMPILER)
add_executable(transmag ${LOOP_SRC})
target_link_libraries(transmag alps Boost::mpi ${ALPS_EXTRA_LIBRARIES})

add_alps_test(transmag)
configure_file(alternating_chain.xml.in alternating_chain.xml)
