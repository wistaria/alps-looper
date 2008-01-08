dist_noinst_DATA = loop.op
check_PROGRAMS = loop loop_mpi
loop_SOURCES = loop.C observable.h options.h
loop_mpi_SOURCES = loop_mpi.C observable.h options.h
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@

TESTS = loop
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check

archive:
	rm -rf loop-$(PACKAGE_VERSION)
	mkdir loop-$(PACKAGE_VERSION)
	cp -a $(top_srcdir)/LICENSE.txt $(top_srcdir)/CITATIONS.txt $(srcdir)/README.ja $(srcdir)/loop.C $(srcdir)/loop.op $(srcdir)/*.h loop-$(PACKAGE_VERSION)
	tar zcf loop-$(PACKAGE_VERSION).tar.gz loop-$(PACKAGE_VERSION)
	rm -rf loop-$(PACKAGE_VERSION)
