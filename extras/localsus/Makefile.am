DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/extras -DLOOP_CONFIG_HEADER=\"localsus_config.h\"

noinst_PROGRAMS =
check_PROGRAMS =
noinst_PROGRAMS += localsus localsus_evaluate
if HAVE_MPI
noinst_PROGRAMS += localsus_mpi
endif
check_PROGRAMS += localsus_ns

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/path_integral.C

dist_noinst_DATA = localsus_ns.ip localsus_ns.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

dist_localsus_SOURCES = localsus_config.h
nodist_localsus_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_localsus_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       localsus_mpi_CPPFLAGS = @CPPFLAGS_MPI@
       localsus_mpi_LDFLAGS = @LDFLAGS_MPI@
       localsus_mpi_LDADD = @LIBS_MPI@
nodist_localsus_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)
nodist_localsus_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)
