add_definitions(-DLOOP_CONFIG_HEADER="extras/top/top_config.h")

add_library(top_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(top_core ${ALPS_LIBRARIES})

add_executable(top ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(top top_core)

if(LOOPER_BUILD_TESTS)
  add_executable(top_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(top_ns top_core)
  add_test(top_ns
    ${CMAKE_COMMAND}
      -Dcmd=top_ns
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=top_ns.ip
      -Doutput=top_ns.op
      -Dtimeout=600
      -P ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake
    )
endif(LOOPER_BUILD_TESTS)
