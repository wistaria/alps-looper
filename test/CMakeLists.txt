set(TESTS alternating_tensor bond_matrix compatibility edge_type
    filtered_graph find_bridge generate_01 gesvd lattice llsp
    model_parameter parse_vector permutation poisson_distribution
    postorder random_choice sign site_matrix spinmodel_helper syev
    temperature union_find virtual_lattice weight)

foreach(name ${TESTS})
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} ${ALPS_LIBRARIES})
  add_test(${name}
    ${CMAKE_COMMAND}
      -Dcmd=${name}
      -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
      -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
      -Ddllexedir=${PROJECT_BINARY_DIR}/bin
      -Dinput=${name}.ip
      -Doutput=${name}.op
      -Dtimeout=600
      -P ${PROJECT_SOURCE_DIR}/config/CMake/run_test.cmake
  )
endforeach(name)
