add_definitions(-DLOOP_CONFIG_HEADER="extras/corrlen/corrlen_config.h")

add_library(corrlen_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(corrlen_core ${ALPS_LIBRARIES})

add_executable(corrlen ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(corrlen corrlen_core)

if(LOOPER_BUILD_TESTS)
  add_executable(corrlen_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(corrlen_ns corrlen_core)
  add_test(corrlen_ns
    ${CMAKE_COMMAND}
      -Dcmd=corrlen_ns
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=corrlen_ns.ip
      -Doutput=corrlen_ns.op
      -Dtimeout=600
      -P ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake
    )
endif(LOOPER_BUILD_TESTS)
