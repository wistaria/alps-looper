AM_CPPFLAGS = @CPPFLAGS_SGL@ @LAPACK_CPPFLAGS@
AM_LDFLAGS = -L$(top_builddir) @LDFLAGS_SGL@ @LAPACK_LDFLAGS@
LDADD = @LAPACK_LIBS@ @LIBS_SGL@

check_PROGRAMS = amida fill_duration hypercubic index_helper permutation random_choice union_find virtual_graph weight xxz_matrix xxz_model
if HAVE_LAPACK
check_PROGRAMS += syev
endif

dist_noinst_DATA = amida.op fill_duration.ip fill_duration.op hypercubic.ip hypercubic.op index_helper.op permutation.op random_choice.op syev.op union_find.op virtual_graph.op weight.ip weight.op xxz_matrix.ip xxz_matrix.op xxz_model.ip xxz_model.op

TESTS = $(check_PROGRAMS)

TESTS_ENVIRONMENT = $(top_srcdir)/run-test $(srcdir)

amida_SOURCES = amida.C
fill_duration_SOURCES = fill_duration.C
hypercubic_SOURCES = hypercubic.C
index_helper_SOURCES = index_helper.C
permutation_SOURCES = permutation.C
random_choice_SOURCES = random_choice.C
union_find_SOURCES = union_find.C
virtual_graph_SOURCES = virtual_graph.C
weight_SOURCES = weight.C
xxz_matrix_SOURCES = xxz_matrix.C
xxz_model_SOURCES = xxz_model.C
syev_SOURCES = syev.C

clean: clean-am clean-extra
clean-extra:
	rm -f $(chech_PROGRAMS) amida.xdr
