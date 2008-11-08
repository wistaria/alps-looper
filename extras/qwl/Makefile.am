#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir)

#
# programs
#

noinst_PROGRAMS =
#noinst_PROGRAMS += pi_ns sse_ns
noinst_PROGRAMS += sse_ns sse_ps
if HAVE_BOOST_MPI
noinst_PROGRAMS += sse_pp
endif

SRC = $(top_srcdir)/loop_custom.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C

#  dist_pi_ns_SOURCES = pi.C
#nodist_pi_ns_SOURCES = $(SRC)
  dist_sse_ns_SOURCES = sse_ns.C sse.h
nodist_sse_ns_SOURCES = $(SRC)
  dist_sse_ps_SOURCES = sse.C sse.h
nodist_sse_ps_SOURCES = $(SRC)
  dist_sse_pp_SOURCES = sse.C sse.h
nodist_sse_pp_SOURCES = $(SRC)
       sse_pp_LDFLAGS = @LDFLAGS_MPI@
       sse_pp_LDADD = @LIBS_MPI@

#
# tests
#

dist_noinst_DATA = sse_ns.ip sse_ns.op
check_PROGRAMS =
#check_PROGRAMS += pi_ns sse_ns
#check_PROGRAMS += sse_ns

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
