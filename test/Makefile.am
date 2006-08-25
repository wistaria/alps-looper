dist_noinst_DATA = \
	alternating_tensor.op \
	bond_matrix.ip bond_matrix.op \
	compatibility.op \
	edge_type.ip edge_type.op \
	filtered_graph.ip filtered_graph.op \
	find_bridge.ip find_bridge.op \
	gesvd.ip gesvd.op \
	hypercubic.ip hypercubic.op \
	llsp.ip llsp.op \
	model_parameter.ip model_parameter.op \
	parse_vector.ip parse_vector.op \
	permutation.op \
	postorder.ip postorder.op \
	random_choice.op \
	range.ip range.op \
	sign.ip sign.op \
	site_matrix.ip site_matrix.op \
	syev.op \
	temperature.ip temperature.op \
	union_find.op \
	virtual_lattice.op \
	weight.ip weight.op

check_PROGRAMS = \
	alternating_tensor \
	bond_matrix \
	compatibility \
	edge_type \
	filtered_graph \
	find_bridge \
	generate_01 \
	gesvd \
	hypercubic \
	llsp \
	model_parameter \
	parse_vector \
	permutation \
	postorder \
	random_choice \
	range \
	sign \
	site_matrix \
	syev \
	temperature \
	union_find \
	virtual_lattice \
	weight

alternating_tensor_SOURCES = alternating_tensor.C
bond_matrix_SOURCES = bond_matrix.C
compatibility_SOURCES = compatibility.C
edge_type_SOURCES = edge_type.C
filtered_graph_SOURCES = filtered_graph.C
find_bridge_SOURCES = find_bridge.C
generate_01_SOURCES = generate_01.C
gesvd_SOURCES = gesvd.C
hypercubic_SOURCES = hypercubic.C
llsp_SOURCES = llsp.C
model_parameter_SOURCES = model_parameter.C
parse_vector_SOURCES = parse_vector.C
permutation_SOURCES = permutation.C
postorder_SOURCES = postorder.C
random_choice_SOURCES = random_choice.C
range_SOURCES = range.C
sign_SOURCES = sign.C
site_matrix_SOURCES = site_matrix.C
syev_SOURCES = syev.C
temperature_SOURCES = temperature.C
union_find_SOURCES = union_find.C
virtual_lattice_SOURCES = virtual_lattice.C
weight_SOURCES = weight.C

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
