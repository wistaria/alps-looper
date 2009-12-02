#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"corrlen_config.h\"

#
# programs
#

noinst_PROGRAMS = corrlen
if HAVE_MPI
if HAVE_BOOST_MPI
noinst_PROGRAMS += corrlen_mpi
endif
endif

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_corrlen_SOURCES = corrlen_config.h
nodist_corrlen_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_corrlen_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       corrlen_mpi_LDFLAGS = @LDFLAGS_MPI@
       corrlen_mpi_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = corrlen_ns.ip corrlen_ns.op

check_PROGRAMS =
check_PROGRAMS += corrlen_ns

nodist_corrlen_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
