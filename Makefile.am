SUBDIRS = doc test ed .

distdir = alps-looper-$(VERSION)

CPPFLAGS_SGL = @CPPFLAGS_SGL@ @LAPACK_CPPFLAGS@ @BINDINGS_CPPFLAGS@
LDFLAGS_SGL = -L$(top_builddir) @LDFLAGS_SGL@ @LAPACK_LDFLAGS@
LIBS_SGL = @LAPACK_LIBS@ @LIBS_SGL@

CPPFLAGS_MPI = @CPPFLAGS_MPI@ @LAPACK_CPPFLAGS@ @BINDINGS_CPPFLAGS@
LDFLAGS_MPI = -L$(top_builddir) @LDFLAGS_MPI@ @LAPACK_LDFLAGS@
LIBS_MPI = @LAPACK_LIBS@ @LIBS_MPI@

CPPFLAGS_PVM = @CPPFLAGS_PVM@ @LAPACK_CPPFLAGS@ @BINDINGS_CPPFLAGS@
LDFLAGS_PVM = -L$(top_builddir) @LDFLAGS_PVM@ @LAPACK_LDFLAGS@
LIBS_PVM = @LAPACK_LIBS@ @LIBS_PVM@

noinst_PROGRAMS = loop evaluate loop_cmd
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
if HAVE_PVM
noinst_PROGRAMS += loop_pvm
endif

nodist_pkginclude_HEADERS = looper/config.h
pkginclude_HEADERS = looper/amida.h looper/copyright.h looper/exact_diag.h looper/fill_duration.h looper/graph.h looper/lapack.h looper/measurement.h looper/model.h looper/node.h looper/operator.h looper/path_integral.h looper/permutation.h looper/random_choice.h looper/sse.h looper/union_find.h looper/util.h looper/vector_helper.h looper/virtual_graph.h looper/weight.h looper/xxz.h

dist_noinst_SCRIPTS = config/run-test config/update_preamble
dist_noinst_DATA = CITATIONS.txt LICENSE lattices.xml models.xml

loop_SOURCES = loop.C loop_impl.h
loop_CPPFLAGS = $(CPPFLAGS_SGL)
loop_LDFLAGS = $(LDFLAGS_SGL)
loop_LDADD = $(LIBS_SGL)

evaluate_SOURCES = evaluate.C loop_impl.h
evaluate_CPPFLAGS = $(CPPFLAGS_SGL)
evaluate_LDFLAGS = $(LDFLAGS_SGL)
evaluate_LDADD = $(LIBS_SGL)

loop_mpi_SOURCES = loop.C loop_impl.h
loop_mpi_CPPFLAGS = $(CPPFLAGS_MPI)
loop_mpi_LDFLAGS = $(LDFLAGS_MPI)
loop_mpi_LDADD = $(LIBS_MPI)

loop_pvm_SOURCES = loop.C loop_impl.h
loop_pvm_CPPFLAGS = $(CPPFLAGS_PVM)
loop_pvm_LDFLAGS = $(LDFLAGS_PVM)
loop_pvm_LDADD = $(LIBS_PVM)

loop_cmd_SOURCES = loop_cmd.C loop_impl.h
loop_cmd_CPPFLAGS = $(CPPFLAGS_SGL)
loop_cmd_LDFLAGS = $(LDFLAGS_SGL)
loop_cmd_LDADD = $(LIBS_SGL)

update-preamble:
	if test -f config/update_preamble; then \
	  find . -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs ./config/update_preamble; \
	fi
