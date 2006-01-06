DEFAULT_INCLUDES = -I$(top_builddir)/v3.1 -I$(top_srcdir)/v3.1

PROG = loop loop_evaluate loop_cmd

noinst_PROGRAMS =
if BUILD_V31
noinst_PROGRAMS += $(PROG)
endif

dist_noinst_HEADERS = looper.h \
	looper/amida.h \
	looper/exact_diag.h \
	looper/fill_duration.h \
	looper/find_bridges.h \
	looper/graph.h \
	looper/lapack.h \
	looper/measurement.h \
	looper/model.h \
	looper/node.h \
	looper/operator.h \
	looper/path_integral.h \
	looper/permutation.h \
	looper/random_choice.h \
	looper/sign.h \
	looper/sse.h \
	looper/union_find.h \
	looper/util.h \
	looper/vector_helper.h \
	looper/weight.h

loop_SOURCES = loop.C loop_factory.C loop_factory.h loop_worker.h
loop_evaluate_SOURCES = loop_evaluate.C loop_worker.h
loop_cmd_SOURCES = loop_cmd.C loop_worker.h

prog: $(PROG)
clean-local:
	rm -f $(PROG)
