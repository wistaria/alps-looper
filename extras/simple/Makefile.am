DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"simple_config.h\"

# programs
noinst_PROGRAMS = loop
  dist_loop_SOURCES = simple_config.h
nodist_loop_SOURCES = $(top_srcdir)/loop.C $(top_srcdir)/path_integral.C $(top_srcdir)/sse.C

# tests
dist_noinst_DATA = loop_ns.ip loop_ns.op
check_PROGRAMS = loop_ns
nodist_loop_ns_SOURCES = $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
