#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir)

#
# programs
#

noinst_PROGRAMS =

#
# tests
#

dist_noinst_DATA =
check_PROGRAMS =

dist_noinst_DATA += ising_ns.ip ising_ns.op
check_PROGRAMS += ising_ns
nodist_ising_ns_SOURCES = $(top_srcdir)/loop_ns.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/ising.C ising_config.h
       ising_ns_CPPFLAGS = @CPPFLAGS@ -DLOOP_CONFIG_HEADER=\"ising_config.h\"

dist_noinst_DATA += xyz_ns.ip xyz_ns.op
check_PROGRAMS += xyz_ns
nodist_xyz_ns_SOURCES = $(top_srcdir)/loop_ns.C $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/diag.C $(top_srcdir)/path_integral.C $(top_srcdir)/sse.C xyz_config.h
       xyz_ns_CPPFLAGS = @CPPFLAGS@ -DLOOP_CONFIG_HEADER=\"xyz_config.h\"

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
