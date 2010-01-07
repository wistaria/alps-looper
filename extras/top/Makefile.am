#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"top_config.h\"

#
# programs
#

noinst_PROGRAMS = top

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_top_SOURCES = top_config.h
nodist_top_SOURCES = $(top_srcdir)/loop.C $(SRC)

#
# tests
#

dist_noinst_DATA = top_ns.ip top_ns.op

check_PROGRAMS =
check_PROGRAMS += top_ns

nodist_top_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
