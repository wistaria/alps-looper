#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"xyz_config.h\"

#
# programs
#

noinst_PROGRAMS = xyz xyz_evaluate
if HAVE_MPI
noinst_PROGRAMS += xyz_mpi
endif
if HAVE_PARAPACK
noinst_PROGRAMS += xyz_ps xyz_pevaluate
if HAVE_MPI
noinst_PROGRAMS += xyz_pp
endif
endif

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/diag.C $(top_srcdir)/path_integral.C $(top_srcdir)/sse.C

  dist_xyz_SOURCES = xyz_config.h
nodist_xyz_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_xyz_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       xyz_mpi_LDFLAGS = @LDFLAGS_MPI@
       xyz_mpi_LDADD = @LIBS_MPI@
nodist_xyz_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)

nodist_xyz_ps_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
nodist_xyz_pp_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
       xyz_pp_LDFLAGS = @LDFLAGS_MPI@
       xyz_pp_LDADD = @LIBS_MPI@
nodist_xyz_pevaluate_SOURCES = $(top_srcdir)/loop_pevaluate.C $(SRC)

#
# tests
#

dist_noinst_DATA = xyz_ns.ip xyz_ns.op

check_PROGRAMS =
if HAVE_PARAPACK
check_PROGRAMS += xyz_ns
endif

nodist_xyz_ns_SOURCES = $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
