set(TESTS alternating_tensor bond_matrix compatibility edge_type
    filtered_graph find_bridge generate_01 gesvd lattice llsp
    model_parameter parse_vector permutation poisson_distribution
    postorder random_choice sign site_matrix spinmodel_helper syev
    temperature union_find virtual_lattice weight)

foreach(name ${TESTS})
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} ${ALPS_LIBRARIES})
  add_alps_test(${name})
endforeach(name)
