SUBDIRS = doc
if EXIST_CHECK
SUBDIRS += check
endif
if EXIST_TEST
SUBDIRS += test
endif
if EXIST_V31
SUBDIRS += v3.1
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

REMOVE_FOR_RELEASE = check

distdir = alps-looper-$(VERSION)

noinst_PROGRAMS = loop loop_evaluate
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
check_PROGRAMS = diag loop_ns

include_HEADERS = looper.h
pkginclude_HEADERS = \
	looper/cluster.h \
	looper/evaluate.h \
	looper/find_bridge.h \
	looper/graph.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/measurement.h \
	looper/model.h \
	looper/montecarlo.h \
	looper/operator.h \
	looper/permutation.h \
	looper/qwl_histogram.h \
	looper/random_choice.h \
	looper/time.h \
	looper/type.h \
	looper/union_find.h \
	looper/util.h \
	looper/version.h \
	looper/weight.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE.txt diag.ip diag.op loop_ns.ip loop_ns.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

SRC = loop_factory.h loop_factory.C \
	loop_worker.h loop_worker.C \
	loop_worker_pi.C \
	loop_worker_sse.C

loop_SOURCES = loop.C loop_config.h $(SRC)
loop_evaluate_SOURCES = loop_evaluate.C loop_config.h
loop_mpi_SOURCES = loop.C loop_config.h $(SRC)
loop_mpi_CPPFLAGS = @CPPFLAGS_MPI@
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

loop_ns_SOURCES = loop_ns.C loop_config.h $(SRC)
diag_SOURCES = diag.C

.PHONY : test tests
test tests : check

.PHONY : update-preamble
update-preamble:
	if test -x "$(top_srcdir)/config/update_preamble"; then \
	  find $(top_srcdir) -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs $(top_srcdir)/config/update_preamble; \
	fi

.PHONY : release snapshot
release snapshot :
	rm -rf $(distdir)
	$(MAKE) dist
	tar zxf $(distdir).tar.gz
	rm -f $(distdir).tar.gz
	cd $(distdir); rm -rf $(REMOVE_FOR_RELEASE)
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
	  rm -rf alps-looper-$$date; \
	  mv -f $(distdir) alps-looper-$$date; \
	  tar zcf alps-looper-$$date.tar.gz alps-looper-$$date; \
	  rm -rf alps-looper-$$date; \
	else \
	  tar zcf $(distdir).tar.gz $(distdir); \
	  rm -rf $(distdir); \
	fi
