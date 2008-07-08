#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"parallel_config.h\"

#
# programs
#

noinst_PROGRAMS =
if HAVE_PARAPACK
if HAVE_MPI
noinst_PROGRAMS += parallel_pm parallel_pevaluate
endif
endif

dist_parallel_pm_SOURCES = loop_pm.C parallel_config.h loop_factory_mpi.h loop_factory_mpi.C path_integral_mpi.C
nodist_parallel_pm_SOURCES = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C 
       parallel_pm_LDFLAGS = @LDFLAGS_MPI@
       parallel_pm_LDADD = @LIBS_MPI@
nodist_parallel_pevaluate_SOURCES = $(top_srcdir)/loop_pevaluate.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

#
# tests
#

dist_noinst_DATA = parallel_ns.ip parallel_ns.op

check_PROGRAMS =
if HAVE_MPI
check_PROGRAMS += parallel_ns
endif

dist_parallel_ns_SOURCES = parallel_ns.C parallel_config.h loop_factory_mpi.h loop_factory_mpi.C path_integral_mpi.C
nodist_parallel_ns_SOURCES = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C 
       parallel_ns_LDFLAGS = @LDFLAGS_MPI@
       parallel_ns_LDADD = @LIBS_MPI@

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
