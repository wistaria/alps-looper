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

include_HEADERS = looper.h
pkginclude_HEADERS = \
	looper/amida.h \
	looper/exact_diag.h \
	looper/fill_duration.h \
	looper/find_bridges.h \
	looper/lapack.h \
	looper/lattice.h \
	looper/measurement.h \
	looper/model.h \
	looper/node.h \
	looper/operator.h \
	looper/path_integral.h \
	looper/permutation.h \
	looper/qmc.h \
	looper/random_choice.h \
	looper/sign.h \
	looper/sse.h \
	looper/type.h \
	looper/union_find.h \
	looper/util.h \
	looper/vector_helper.h \
	looper/weight.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE.txt

loop_SOURCES = loop.C loop_factory.C loop_factory.h loop_worker.h
loop_evaluate_SOURCES = loop_evaluate.C loop_worker.h
loop_nosched_SOURCES = loop.C loop_factory.C loop_factory.h loop_worker.h
loop_nosched_CPPFLAGS = -DWITHOUT_SCHEDULER

loop_mpi_SOURCES = loop.C loop_factory.C loop_factory.h loop_worker.h
loop_mpi_CPPFLAGS = @CPPFLAGS_MPI@
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

loop_pvm_SOURCES = loop.C loop_factory.C loop_factory.h loop_worker.h
loop_pvm_CPPFLAGS = @CPPFLAGS_PVM@
loop_pvm_LDFLAGS = @LDFLAGS_PVM@
loop_pvm_LDADD = @LIBS_PVM@

.PHONY : test tests
test tests :
	$(MAKE) -C test $@

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
