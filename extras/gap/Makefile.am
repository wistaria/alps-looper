DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/extras

noinst_PROGRAMS =
check_PROGRAMS =
noinst_PROGRAMS += gap gap_evaluate
if HAVE_MPI
noinst_PROGRAMS += gap_mpi
endif
check_PROGRAMS += gap_ns

SRC = $(top_srcdir)/loop.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_worker.C $(top_srcdir)/loop_worker_pi.C
SRC_EVALUATE = $(top_srcdir)/loop_evaluate.C
SRC_NS = $(top_srcdir)/loop_ns.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_worker.C $(top_srcdir)/loop_worker_pi.C

dist_noinst_DATA = gap_ns.ip gap_ns.op

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

dist_gap_SOURCES = gap_config.h
nodist_gap_SOURCES = $(SRC)
       gap_CPPFLAGS = @CPPFLAGS@ -DLOOP_CONFIG_HEADER=\"gap_config.h\"
nodist_gap_mpi_SOURCES = $(SRC)
       gap_mpi_CPPFLAGS = @CPPFLAGS_MPI@ -DLOOP_CONFIG_HEADER=\"gap_config.h\"
       gap_mpi_LDFLAGS = @LDFLAGS_MPI@
       gap_mpi_LDADD = @LIBS_MPI@
nodist_gap_evaluate_SOURCES = $(SRC_EVALUTE)
       gap_evaluate_CPPFLAGS = @CPPFLAGS@ -DLOOP_CONFIG_HEADER=\"gap_config.h\"
nodist_gap_ns_SOURCES = $(SRC_NS)
       gap_ns_CPPFLAGS = @CPPFLAGS@ -DLOOP_CONFIG_HEADER=\"gap_config.h\"
