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
	looper/crop.h \
	looper/divide_if_positive.h \
	looper/evaluator.h \
	looper/evaluator_impl.h \
	looper/find_bridge.h \
	looper/flatten_matrix.h \
	looper/gap.h \
	looper/generate_seed.h \
	looper/graph.h \
	looper/histogram.h \
	looper/integer_range.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/localsus.h \
	looper/measurement.h \
	looper/model.h \
	looper/model_impl.h \
	looper/montecarlo.h \
	looper/operator.h \
	looper/permutation.h \
	looper/power.h \
	looper/random_choice.h \
	looper/stiffness.h \
	looper/susceptibility.h \
	looper/temperature.h \
	looper/top.h \
	looper/time.h \
	looper/transmag.h \
	looper/type.h \
	looper/union_find.h \
	looper/version.h \
	looper/weight.h

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
noinst_PROGRAMS += loop_pp
endif

SRC = loop_config.h loop_factory.h loop_factory.C loop_model.C
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

loop_evaluate_SOURCES = loop_evaluate.C qwl_evaluate.C $(SRC)

loop_pp_SOURCES = loop_pp.C $(SRC)
loop_pp_LDFLAGS = @LDFLAGS_PP@
loop_pp_LDADD = @LIBS_PP@

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

REMOVE_FOR_STANDARD = check extras looper/extras test
REMOVE_FOR_PRO = check

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
	  rm -rf $(distdir); \
	  mv -f $(distdir) $(distdir); \
	  tar zcf $(distdir).tar.gz $(distdir); \
	  mv -f $(distdir) $(distdir); \
	fi
	rm -rf $(distdir)
