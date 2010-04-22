add_definitions(-DLOOP_CONFIG_HEADER="extras/top/top_config.h")

add_library(top_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(top_core ${ALPS_LIBRARIES})

add_executable(top ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(top top_core)

if(LOOPER_BUILD_TESTS)
  add_executable(top_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(top_ns top_core)
  add_alps_test(top_ns)
endif(LOOPER_BUILD_TESTS)
