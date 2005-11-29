SUBDIRS = doc test ed check .

REMOVE_FOR_RELEASE =

distdir = alps-looper-$(VERSION)

noinst_PROGRAMS = loop loop_evaluate loop_nosched
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
if HAVE_PVM
noinst_PROGRAMS += loop_pvm
endif
check_PROGRAMS = loop_nosched

include_HEADERS = looper.h
pkginclude_HEADERS = \
	looper/find_bridges.h \
	looper/graph.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/location.h \
	looper/model.h \
	looper/operator.h \
	looper/permutation.h \
	looper/random_choice.h \
	looper/type.h \
	looper/union_find.h \
	looper/util.h \
	looper/weight.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE.txt loop_nosched.ip loop_nosched.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

SRC = loop.C loop_config.h loop_factory.C loop_factory.h \
	loop_worker.h loop_worker.C \
	loop_worker_pi.h loop_worker_pi.C \
	loop_worker_sse.h loop_worker_sse.C

loop_SOURCES = $(SRC)
loop_evaluate_SOURCES = loop_evaluate.C
loop_nosched_SOURCES = $(SRC)
loop_nosched_CPPFLAGS = -DWITHOUT_SCHEDULER

loop_mpi_SOURCES = $(SRC)
loop_mpi_CPPFLAGS = @CPPFLAGS_MPI@
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

loop_pvm_SOURCES = $(SRC)
loop_pvm_CPPFLAGS = @CPPFLAGS_PVM@
loop_pvm_LDFLAGS = @LDFLAGS_PVM@
loop_pvm_LDADD = @LIBS_PVM@

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
