check_PROGRAMS = \
	alternating_tensor \
	amida \
	bond_matrix \
	edge_type \
	fill_duration \
	filtered_graph \
	find_bridges \
	generate_01 \
	gesvd \
	hypercubic \
	index_helper \
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
	union_find \
	virtual_graph \
	weight

dist_noinst_DATA = \
	alternating_tensor.op \
	amida.op \
	bond_matrix.ip bond_matrix.op \
	edge_type.ip edge_type.op \
	fill_duration.ip fill_duration.op \
	filtered_graph.ip filtered_graph.op \
	find_bridges.ip find_bridges.op \
	gesvd.ip gesvd.op \
	hypercubic.ip hypercubic.op \
	index_helper.op \
	llsp.ip llsp.op \
	model_parameter.op \
	parse_vector.ip parse_vector.op \
	permutation.op \
	postorder.ip postorder.op \
	random_choice.op \
	range.ip range.op \
	sign.ip sign.op \
	site_matrix.ip site_matrix.op \
	syev.op \
	union_find.op \
	virtual_graph.op \
	weight.ip weight.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

alternating_tensor_SOURCES = alternating_tensor.C
amida_SOURCES = amida.C
bond_matrix_SOURCES = bond_matrix.C
edge_type_SOURCES = edge_type.C
fill_duration_SOURCES = fill_duration.C
filtered_graph_SOURCES = filtered_graph.C
find_bridges_SOURCES = find_bridges.C
generate_01_SOURCES = generate_01.C
gesvd_SOURCES = gesvd.C
hypercubic_SOURCES = hypercubic.C
index_helper_SOURCES = index_helper.C
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
union_find_SOURCES = union_find.C
virtual_graph_SOURCES = virtual_graph.C
weight_SOURCES = weight.C

test tests : check
clean-local :
	rm -f amida.xdr
