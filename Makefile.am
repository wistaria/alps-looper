SUBDIRS = doc test ed

CPPFLAGS_SGL = @CPPFLAGS_SGL@ @LAPACK_CPPFLAGS@
LDFLAGS_SGL = -L$(top_builddir) @LDFLAGS_SGL@ @LAPACK_LDFLAGS@
LIBS_SGL = @LAPACK_LIBS@ @LIBS_SGL@

CPPFLAGS_MPI = @CPPFLAGS_MPI@ @LAPACK_CPPFLAGS@
LDFLAGS_MPI = -L$(top_builddir) @LDFLAGS_MPI@ @LAPACK_LDFLAGS@
LIBS_MPI = @LAPACK_LIBS@ @LIBS_MPI@

CPPFLAGS_PVM = @CPPFLAGS_PVM@ @LAPACK_CPPFLAGS@
LDFLAGS_PVM = -L$(top_builddir) @LDFLAGS_PVM@ @LAPACK_LDFLAGS@
LIBS_PVM = @LAPACK_LIBS@ @LIBS_PVM@

bin_PROGRAMS = qmc qmc_cmd
if HAVE_MPI
bin_PROGRAMS += qmc_mpi
endif
if HAVE_PVM
bin_PROGRAMS += qmc_pvm
endif

nodist_pkginclude_HEADERS = looper/config.h
pkginclude_HEADERS = looper/amida.h looper/copyright.h looper/exact_diag.h looper/fill_duration.h looper/fortran.h looper/graph.h looper/lapack.h looper/measurement.h looper/model.h looper/node.h looper/path_integral.h looper/permutation.h looper/random_choice.h looper/sse.h looper/union_find.h looper/util.h looper/vector_helper.h looper/virtual_graph.h looper/weight.h looper/xxz.h

dist_noinst_SCRIPTS = run-test

qmc_SOURCES = qmc.C qmc_impl.h
qmc_CPPFLAGS = $(CPPFLAGS_SGL)
qmc_LDFLAGS = $(LDFLAGS_SGL)
qmc_LDADD = $(LIBS_SGL)

qmc_mpi_SOURCES = qmc.C qmc_impl.h
qmc_mpi_CPPFLAGS = $(CPPFLAGS_MPI)
qmc_mpi_LDFLAGS = $(LDFLAGS_MPI)
qmc_mpi_LDADD = $(LIBS_MPI)

qmc_pvm_SOURCES = qmc.C qmc_impl.h
qmc_pvm_CPPFLAGS = $(CPPFLAGS_PVM)
qmc_pvm_LDFLAGS = $(LDFLAGS_PVM)
qmc_pvm_LDADD = $(LIBS_PVM)

qmc_cmd_SOURCES = qmc_cmd.C qmc_impl.h
qmc_cmd_CPPFLAGS = $(CPPFLAGS_SGL)
qmc_cmd_LDFLAGS = $(LDFLAGS_SGL)
qmc_cmd_LDADD = $(LIBS_SGL)

update-preamble:
	find . -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs ./update_preamble.pl
