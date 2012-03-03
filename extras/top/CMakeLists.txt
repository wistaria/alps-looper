add_definitions(-DLOOP_CONFIG_HEADER="extras/top/top_config.h")

add_executable(top ${PROJECT_SOURCE_DIR}/loop.C ${PROJECT_SOURCE_DIR}/loop_model.C
  ${PROJECT_SOURCE_DIR}/path_integral.C)
target_link_libraries(top ${ALPS_LIBRARIES})

add_alps_test(top)
