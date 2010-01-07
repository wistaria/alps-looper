#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"localsus_config.h\"

#
# programs
#

noinst_PROGRAMS = localsus

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_localsus_SOURCES = localsus_config.h
nodist_localsus_SOURCES = $(top_srcdir)/loop.C $(SRC)

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
