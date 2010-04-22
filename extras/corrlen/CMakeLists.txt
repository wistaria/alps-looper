add_definitions(-DLOOP_CONFIG_HEADER="extras/corrlen/corrlen_config.h")

add_library(corrlen_core ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(corrlen_core ${ALPS_LIBRARIES})

add_executable(corrlen ${PROJECT_SOURCE_DIR}/loop.C)
target_link_libraries(corrlen corrlen_core)

if(LOOPER_BUILD_TESTS)
  add_executable(corrlen_ns ${PROJECT_SOURCE_DIR}/loop_ns.C)
  target_link_libraries(corrlen_ns corrlen_core)
  add_alps_test(corrlen_ns)
endif(LOOPER_BUILD_TESTS)
