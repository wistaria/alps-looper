if(MPI_FOUND)
  add_definitions(-DLOOP_CONFIG_HEADER="extras/parallel/parallel_config.h")

  add_library(parallel_core ${PROJECT_SOURCE_DIR}/loop_model.C
    ${PROJECT_SOURCE_DIR}/path_integral.C path_integral_mpi.C)
  target_link_libraries(parallel_core ${ALPS_LIBRARIES})

  add_executable(parallel ${PROJECT_SOURCE_DIR}/loop.C)
  target_link_libraries(parallel parallel_core)

  add_executable(parallel_ns parallel_ns.C)
  target_link_libraries(parallel_ns parallel_core)

  configure_file(test.sh.in test.sh)
endif(MPI_FOUND)
