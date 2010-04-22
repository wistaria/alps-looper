add_definitions(-DLOOP_CONFIG_HEADER="extras/xyz/xyz_config.h")

add_library(xyz_core ${PROJECT_SOURCE_DIR}/loop_model.C ${PROJECT_SOURCE_DIR}/diag.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(xyz_core ${ALPS_LIBRARIES})

add_executable(xyz ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(xyz xyz_core)

if(LOOPER_BUILD_TESTS)
  add_executable(xyz_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(xyz_ns xyz_core)
  add_test(xyz_ns
    ${CMAKE_COMMAND}
      -Dcmd=xyz_ns
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=xyz_ns.ip
      -Doutput=xyz_ns.op
      -Dtimeout=600
      -P ${CMAKE_INSTALL_PREFIX}/share/alps/run_test.cmake
    )
  configure_file(xyz_model.xml.in xyz_model.xml)
endif(LOOPER_BUILD_TESTS)
