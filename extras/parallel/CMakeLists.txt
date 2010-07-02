if(MPI_FOUND)
  add_definitions(-DLOOP_CONFIG_HEADER="extras/parallel/parallel_config.h")

  add_executable(parallel ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
    ${PROJECT_SOURCE_DIR}/path_integral.C path_integral_mpi.C)
  target_link_libraries(parallel ${ALPS_LIBRARIES})

  configure_file(test.sh.in test.sh)
endif(MPI_FOUND)
