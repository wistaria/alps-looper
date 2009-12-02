#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"xyz_config.h\"

#
# programs
#

noinst_PROGRAMS = xyz
if HAVE_MPI
if HAVE_BOOST_MPI
noinst_PROGRAMS += xyz_mpi
endif
endif

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/diag.C $(top_srcdir)/path_integral.C $(top_srcdir)/sse.C

  dist_xyz_SOURCES = xyz_config.h
nodist_xyz_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_xyz_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       xyz_mpi_LDFLAGS = @LDFLAGS_MPI@
       xyz_mpi_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = xyz_ns.ip xyz_ns.op

check_PROGRAMS =
check_PROGRAMS += xyz_ns

nodist_xyz_ns_SOURCES = $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
