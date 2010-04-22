add_definitions(-DLOOP_CONFIG_HEADER="extras/gap/gap_config.h")

add_library(gap_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(gap_core ${ALPS_LIBRARIES})

add_executable(gap ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(gap gap_core)

if(LOOPER_BUILD_TESTS)
  add_executable(gap_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(gap_ns gap_core)
  add_test(gap_ns
    ${CMAKE_COMMAND}
      -Dcmd=gap_ns
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=gap_ns.ip
      -Doutput=gap_ns.op
      -Dtimeout=600
      -P ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake
    )
endif(LOOPER_BUILD_TESTS)
