#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir)

#
# programs
#

noinst_PROGRAMS =
if HAVE_PARAPACK
#noinst_PROGRAMS += pi_ns sse_ns
noinst_PROGRAMS += sse_ns
if HAVE_MPI
noinst_PROGRAMS +=
endif
endif

SRC = $(top_srcdir)/loop_custom.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C

#  dist_pi_ns_SOURCES = pi.C
#nodist_pi_ns_SOURCES = $(SRC)
  dist_sse_ns_SOURCES = sse_ns.C sse.h
nodist_sse_ns_SOURCES = $(SRC)

#
# tests
#

dist_noinst_DATA =
check_PROGRAMS =
if HAVE_PARAPACK
#check_PROGRAMS += pi_ns sse_ns
check_PROGRAMS += sse_ns
endif

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
