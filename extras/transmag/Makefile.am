DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/extras -DLOOP_CONFIG_HEADER=\"transmag_config.h\"

noinst_PROGRAMS =
check_PROGRAMS =
noinst_PROGRAMS += transmag transmag_evaluate
if HAVE_MPI
noinst_PROGRAMS += transmag_mpi
endif
check_PROGRAMS += transmag_ns

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/diag.C $(top_srcdir)/path_integral.C

dist_noinst_DATA = transmag_ns.ip transmag_ns.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

dist_transmag_SOURCES = transmag_config.h
nodist_transmag_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_transmag_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       transmag_mpi_CPPFLAGS = @CPPFLAGS_MPI@
       transmag_mpi_LDFLAGS = @LDFLAGS_MPI@
       transmag_mpi_LDADD = @LIBS_MPI@
nodist_transmag_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)
nodist_transmag_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)
