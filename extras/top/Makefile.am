#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"top_config.h\"

#
# programs
#

noinst_PROGRAMS = top top_evaluate
if HAVE_MPI
noinst_PROGRAMS += top_mpi
endif
if HAVE_PARAPACK
noinst_PROGRAMS += top_ps top_pevaluate
if HAVE_MPI
noinst_PROGRAMS += top_pp
endif
endif

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_top_SOURCES = top_config.h
nodist_top_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_top_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       top_mpi_LDFLAGS = @LDFLAGS_MPI@
       top_mpi_LDADD = @LIBS_MPI@
nodist_top_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)

nodist_top_ps_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
nodist_top_pp_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
       top_pp_LDFLAGS = @LDFLAGS_MPI@
       top_pp_LDADD = @LIBS_MPI@
nodist_top_pevaluate_SOURCES = $(top_srcdir)/loop_pevaluate.C $(SRC)

#
# tests
#

dist_noinst_DATA = top_ns.ip top_ns.op

check_PROGRAMS = top_ns

nodist_top_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
