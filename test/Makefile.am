check_PROGRAMS = \
	alternating_tensor \
	amida \
	edge_type \
	fill_duration \
	filtered_graph \
	find_bridges \
	gels \
	generate_01 \
	hypercubic \
	index_helper \
	permutation \
	postorder \
	random_choice \
	syev \
	union_find \
	virtual_graph \
	weight \
	xxz_matrix \
	xxz_model

dist_noinst_DATA = \
	alternating_tensor.op \
	amida.op \
	edge_type.ip edge_type.op \
	fill_duration.ip fill_duration.op \
	filtered_graph.ip filtered_graph.op \
	find_bridges.ip find_bridges.op \
	gels.ip gels.op \
	generate_01.op \
	hypercubic.ip hypercubic.op \
	index_helper.op \
	permutation.op \
	postorder.ip postorder.op \
	random_choice.op \
	syev.op \
	union_find.op \
	virtual_graph.op \
	weight.ip weight.op \
	xxz_matrix.ip xxz_matrix.op \
	xxz_model.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

alternating_tensor_SOURCES = alternating_tensor.C
amida_SOURCES = amida.C
edge_type_SOURCES = edge_type.C
fill_duration_SOURCES = fill_duration.C
filtered_graph_SOURCES = filtered_graph.C
find_bridges_SOURCES = find_bridges.C
gels_SOURCES = gels.C
generate_01_SOURCES = generate_01.C
hypercubic_SOURCES = hypercubic.C
index_helper_SOURCES = index_helper.C
permutation_SOURCES = permutation.C
postorder_SOURCES = postorder.C
random_choice_SOURCES = random_choice.C
syev_SOURCES = syev.C
union_find_SOURCES = union_find.C
virtual_graph_SOURCES = virtual_graph.C
weight_SOURCES = weight.C
xxz_matrix_SOURCES = xxz_matrix.C
xxz_model_SOURCES = xxz_model.C

test tests : check
clean-local :
	rm -f amida.xdr
