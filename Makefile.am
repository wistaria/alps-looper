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
DIST_SUBDIRS = $(SUBDIRS)

#
# files
#

dist_noinst_HEADERS = \
	looper.h \
	looper/cluster.h \
	looper/evaluate.h \
	looper/find_bridge.h \
	looper/graph.h \
	looper/histogram.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/measurement.h \
	looper/model.h \
	looper/montecarlo.h \
	looper/operator.h \
	looper/permutation.h \
	looper/random_choice.h \
	looper/time.h \
	looper/type.h \
	looper/union_find.h \
	looper/util.h \
	looper/version.h \
	looper/weight.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE.txt diag.ip diag.op loop_ns.ip loop_ns.op

#
# programs
#

noinst_PROGRAMS = loop loop_evaluate
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
check_PROGRAMS = diag loop_ns

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
