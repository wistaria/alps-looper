#
# subdirectories
#

SUBDIRS = doc
if EXIST_CHECK
SUBDIRS += check
endif
if EXIST_EXTRAS
SUBDIRS += extras
endif
if EXIST_STANDALONE
SUBDIRS += standalone
endif
if EXIST_TEST
SUBDIRS += test
endif
if EXIST_V31
SUBDIRS += v3.1
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

#
# files
#

dist_noinst_HEADERS = \
	looper/alternating_tensor.h \
	looper/cluster.h \
	looper/correlation.h \
	looper/custom.h \
	looper/custom_impl.h \
	looper/crop.h \
	looper/divide_if_positive.h \
	looper/evaluator.h \
	looper/evaluator_impl.h \
	looper/find_bridge.h \
	looper/flatten_matrix.h \
	looper/generate_seed.h \
	looper/graph.h \
	looper/graph_impl.h \
	looper/histogram.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/location_impl.h \
	looper/matrix.h \
	looper/measurement.h \
	looper/model.h \
	looper/model_impl.h \
	looper/model_parameter.h \
	looper/montecarlo.h \
	looper/operator.h \
	looper/permutation.h \
	looper/poisson_distribution.h \
	looper/power.h \
	looper/random_choice.h \
	looper/stiffness.h \
	looper/susceptibility.h \
	looper/temperature.h \
	looper/time.h \
	looper/type.h \
	looper/union_find.h \
	looper/version.h \
	looper/weight.h \
	looper/weight_impl.h

HEADERS_PRO = \
	looper/correlation_length.h \
	looper/gap.h \
	looper/localsus.h \
	looper/parallel.h \
	looper/sop.h \
	looper/top.h \
	looper/transmag.h

if LOOPER_PRO
dist_noinst_HEADERS += $(HEADERS_PRO)
endif

dist_noinst_SCRIPTS = bootstrap setup.sh config/extract_date.sh config/extract_version.sh config/run-test config/update_preamble config/update_version
dist_noinst_DATA = CITATIONS.txt LICENSE.txt README config/preamble.in config/version.in

#
# programs
#

noinst_PROGRAMS = loop
if HAVE_PARAPACK_13
if HAVE_BOOST_MPI
noinst_PROGRAMS += loop_mpi
endif
endif

SRC = loop_config.h loop_custom.C loop_model.C
if ENABLE_DIAG
SRC += diag.C
endif
if ENABLE_ISING
SRC += ising.C
endif
if ENABLE_PI
SRC += path_integral.C
endif
if ENABLE_SSE
SRC += sse.C
endif
if ENABLE_SSE_QWL
SRC += sse_qwl.C
endif

loop_SOURCES = loop.C $(SRC)
loop_mpi_SOURCES = loop.C $(SRC)
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA += loop_ns.ip loop_ns.op

check_PROGRAMS =
if ENABLE_DIAG
if ENABLE_ISING
if ENABLE_PI
if ENABLE_SSE
check_PROGRAMS += loop_ns
endif
endif
endif
endif

loop_ns_SOURCES = loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check

check-extras :
	test -f extras/Makefile && make -C extras check-extras

#
# maintenance
#

.PHONY : update-preamble
update-preamble:
	if test -x "$(top_srcdir)/config/update_preamble"; then \
	  find $(top_srcdir) -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs $(top_srcdir)/config/update_preamble; \
	fi

#
# release and snapshot
#

REMOVE_FOR_STANDARD = check extras looper/extras standalone test sse_qwl.C qwl_evaluate.C v3.1 $(HEADERS_PRO)
REMOVE_FOR_PRO =

distdir = alps-looper-$(VERSION)
distdir_pro = alps-looper-pro-$(VERSION)

.PHONY : release snapshot
release snapshot :
	rm -rf $(distdir)
	$(MAKE) dist
	tar zxf $(distdir).tar.gz
	rm -f $(distdir).tar.gz
	cd $(distdir); rm -rf $(REMOVE_FOR_PRO)
	rm -rf $(distdir_pro)
	mv -f $(distdir) $(distdir_pro)
	tar zcf $(distdir_pro).tar.gz $(distdir_pro)
	mv -f $(distdir_pro) $(distdir)
	cd $(distdir); rm -rf $(REMOVE_FOR_STANDARD)
	tar zcf $(distdir).tar.gz $(distdir)
	rm -rf $(distdir)
