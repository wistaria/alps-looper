SUBDIRS = doc test ed .

distdir = alps-looper-$(VERSION)

noinst_PROGRAMS = loop evaluate loop_cmd
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
if HAVE_PVM
noinst_PROGRAMS += loop_pvm
endif

nodist_pkginclude_HEADERS = looper/config.h
pkginclude_HEADERS = \
	looper/amida.h \
	looper/copyright.h \
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
	looper/sse.h \
	looper/union_find.h \
	looper/util.h \
	looper/vector_helper.h \
	looper/version.h \
	looper/virtual_graph.h \
	looper/weight.h \
	looper/xxz.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE lattices.xml models.xml

loop_SOURCES = loop.C loop_impl.h
evaluate_SOURCES = evaluate.C
loop_cmd_SOURCES = loop_cmd.C

loop_mpi_SOURCES = loop.C
loop_mpi_CPPFLAGS = @CPPFLAGS_MPI@
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

loop_pvm_SOURCES = loop.C
loop_pvm_CPPFLAGS = @CPPFLAGS_PVM@
loop_pvm_LDFLAGS = @LDFLAGS_PVM@
loop_pvm_LDADD = @LIBS_PVM@

.PHONY : test tests
test tests :
	$(MAKE) -C test check

.PHONY : update-preamble
update-preamble:
	if test -f config/update_preamble; then \
	  find . -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs ./config/update_preamble; \
	fi

.PHONY : release snapshot
release : dist
snapshot :
	$(MAKE) dist
	tar zxf $(distdir).tar.gz
	date=`date +%Y%m%d`; \
        mv -f $(distdir) alps-looper-$$date; \
	tar zcf alps-looper-$$date.tar.gz alps-looper-$$date; \
	rm -rf alps-looper-$$date
