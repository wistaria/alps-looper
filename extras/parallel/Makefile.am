#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"parallel_config.h\"

#
# programs
#

noinst_PROGRAMS =
if HAVE_BOOST_MPI
noinst_PROGRAMS += parallel_mpi
endif

dist_parallel_mpi_SOURCES = $(top_srcdir)/loop.C parallel_config.h $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C path_integral_mpi.C

#
# tests
#

dist_noinst_DATA = parallel_ns.ip parallel_ns.op-2 parallel_ns.op-3 parallel_ns.op-4

check_PROGRAMS =
if HAVE_BOOST_MPI
check_PROGRAMS += parallel_ns
endif

dist_parallel_ns_SOURCES = parallel_ns.C parallel_config.h $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C path_integral_mpi.C

TESTS =
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
