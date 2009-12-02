#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"gap_config.h\"

#
# programs
#

noinst_PROGRAMS = gap
if HAVE_MPI
if HAVE_BOOST_MPI
noinst_PROGRAMS += gap_mpi
endif
endif

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_gap_SOURCES = gap_config.h
nodist_gap_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_gap_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       gap_mpi_LDFLAGS = @LDFLAGS_MPI@
       gap_mpi_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = gap_ns.ip gap_ns.op

check_PROGRAMS =
check_PROGRAMS += gap_ns

nodist_gap_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
