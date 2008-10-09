#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"localsus_config.h\"

#
# programs
#

noinst_PROGRAMS = localsus localsus_evaluate
if HAVE_MPI
noinst_PROGRAMS += localsus_mpi
endif
if HAVE_PARAPACK
noinst_PROGRAMS += localsus_ps localsus_pevaluate
if HAVE_MPI
noinst_PROGRAMS += localsus_pp
endif
endif

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_localsus_SOURCES = localsus_config.h
nodist_localsus_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_localsus_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       localsus_mpi_LDFLAGS = @LDFLAGS_MPI@
       localsus_mpi_LDADD = @LIBS_MPI@
nodist_localsus_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)

nodist_localsus_ps_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
nodist_localsus_pp_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
       localsus_pp_LDFLAGS = @LDFLAGS_MPI@
       localsus_pp_LDADD = @LIBS_MPI@
nodist_localsus_pevaluate_SOURCES = $(top_srcdir)/loop_pevaluate.C $(SRC)

#
# tests
#

dist_noinst_DATA = localsus_ns.ip localsus_ns.op

check_PROGRAMS =
if HAVE_PARAPACK
check_PROGRAMS += localsus_ns
endif

nodist_localsus_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
