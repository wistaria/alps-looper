#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"parallel_config.h\"

#
# programs
#

noinst_PROGRAMS = parallel_ps
if HAVE_BOOST_MPI
noinst_PROGRAMS += parallel_pm
endif

dist_parallel_ps_SOURCES = $(top_srcdir)/loop_pp.C parallel_config.h $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

dist_parallel_pm_SOURCES = $(top_srcdir)/loop_pp.C parallel_config.h $(top_srcdir)/loop_model.C path_integral_mpi.C
       parallel_pm_LDFLAGS = @LDFLAGS_MPI@
       parallel_pm_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = parallel_ns.ip parallel_ns.op-2 parallel_ns.op-4

check_PROGRAMS =
if HAVE_BOOST_MPI
check_PROGRAMS += parallel_ns
endif

dist_parallel_ns_SOURCES = parallel_ns.C parallel_config.h $(top_srcdir)/loop_model.C path_integral_mpi.C
       parallel_ns_LDFLAGS = @LDFLAGS_MPI@
       parallel_ns_LDADD = @LIBS_MPI@

TESTS =
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
