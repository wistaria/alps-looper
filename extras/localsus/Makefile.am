#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"localsus_config.h\"

#
# programs
#

noinst_PROGRAMS = localsus
if HAVE_PARAPACK_13
if HAVE_BOOST_MPI
noinst_PROGRAMS += localsus_mpi
endif
endif

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_localsus_SOURCES = localsus_config.h
nodist_localsus_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_localsus_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       localsus_mpi_LDFLAGS = @LDFLAGS_MPI@
       localsus_mpi_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = localsus_ns.ip localsus_ns.op

check_PROGRAMS =
check_PROGRAMS += localsus_ns

nodist_localsus_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
