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
	looper/histogram.h \
	looper/integer_range.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
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
	looper/weight.h

HEADERS_PRO = \
	looper/gap.h \
	looper/localsus.h \
	looper/sop.h \
	looper/top.h \
	looper/transmag.h

if LOOPER_PRO
dist_noinst_HEADERS += $(HEADERS_PRO)
endif

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE.txt

#
# programs
#

noinst_PROGRAMS = loop loop_evaluate
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
if HAVE_PARAPACK
noinst_PROGRAMS += loop_ps loop_pevaluate
if HAVE_MPI
noinst_PROGRAMS += loop_pp
endif
endif

SRC = loop_config.h loop_factory.h loop_factory.C loop_custom.C loop_model.C
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

loop_evaluate_SOURCES = loop_evaluate.C $(SRC)
if ENABLE_SSE_QWL
loop_evaluate_SOURCES += qwl_evaluate.C
endif

loop_mpi_SOURCES = loop.C $(SRC)
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

loop_ps_SOURCES = loop_pp.C $(SRC)
loop_ps_LDFLAGS = @LDFLAGS_PS@
loop_ps_LDADD = @LIBS_PS@

loop_pp_SOURCES = loop_pp.C $(SRC)
loop_pp_LDFLAGS = @LDFLAGS_PP@
loop_pp_LDADD = @LIBS_PP@

loop_pevaluate_SOURCES = loop_pevaluate.C $(SRC)
if ENABLE_SSE_QWL
loop_pevaluate_SOURCES += qwl_evaluate.C
endif
loop_pevaluate_LDFLAGS = @LDFLAGS_PS@
loop_pevaluate_LDADD = @LIBS_PS@

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

REMOVE_FOR_STANDARD = check extras looper/extras test sse_qwl.C qwl_evaluate.C $(HEADERS_PRO)
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
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
	  rm -rf alps-looper-pro-$$date; \
	  mv -f $(distdir) alps-looper-pro-$$date; \
	  tar zcf alps-looper-pro-$$date.tar.gz alps-looper-pro-$$date; \
	  mv -f alps-looper-pro-$$date $(distdir); \
	else \
	  rm -rf $(distdir_pro); \
	  mv -f $(distdir) $(distdir_pro); \
	  tar zcf $(distdir_pro).tar.gz $(distdir_pro); \
	  mv -f $(distdir_pro) $(distdir); \
	fi
	cd $(distdir); rm -rf $(REMOVE_FOR_STANDARD)
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
	  rm -rf alps-looper-$$date; \
	  mv -f $(distdir) alps-looper-$$date; \
	  tar zcf alps-looper-$$date.tar.gz alps-looper-$$date; \
	  mv -f alps-looper-$$date $(distdir); \
	else \
	  tar zcf $(distdir).tar.gz $(distdir); \
	fi
	rm -rf $(distdir)
