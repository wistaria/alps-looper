dist_noinst_DATA = loop.op
check_PROGRAMS = loop loop_mpi
loop_SOURCES = loop.C observable.h options.h
loop_mpi_SOURCES = loop_mpi.C observable.h options.h
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

TESTS = loop
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
