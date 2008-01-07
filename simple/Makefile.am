dist_noinst_DATA = loop.op
check_PROGRAMS = loop
loop_SOURCES = loop.C observable.h options.h

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
