add_definitions(-DLOOP_CONFIG_HEADER="extras/localsus/localsus_config.h")

add_library(localsus_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(localsus_core ${ALPS_LIBRARIES})

add_executable(localsus ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(localsus localsus_core)

if(LOOPER_BUILD_TESTS)
  add_executable(localsus_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(localsus_ns localsus_core)
  add_alps_test(localsus_ns)
endif(LOOPER_BUILD_TESTS)
