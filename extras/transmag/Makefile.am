#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"transmag_config.h\"

#
# programs
#

noinst_PROGRAMS = transmag

SRC = $(top_srcdir)/loop_model.C $(top_srcdir)/diag.C $(top_srcdir)/path_integral.C

  dist_transmag_SOURCES = transmag_config.h
nodist_transmag_SOURCES = $(top_srcdir)/loop.C $(SRC)

#
# tests
#

dist_noinst_DATA = transmag_ns.ip transmag_ns.op

check_PROGRAMS =
check_PROGRAMS += transmag_ns

nodist_transmag_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
