DEFAULT_INCLUDES = -I$(srcdir)
dist_noinst_DATA = loop.op loop_mpi.op-1 loop_mpi.op-2 loop_mpi.op-3 loop_mpi.op-4
check_PROGRAMS = loop
if HAVE_PARAPACK_13
else
if HAVE_MPI
check_PROGRAMS += loop_mpi
endif
endif
loop_SOURCES = loop.C common.h observable.h options.h union_find.h
loop_mpi_SOURCES = loop_mpi.C common.h observable.h options.h parallel.h union_find.h

TESTS = loop
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check

mpitest mpicheck : loop_mpi
	mpirun -np 1 loop_mpi | diff $(srcdir)/loop_mpi.op-1 -
	mpirun -np 2 loop_mpi | diff $(srcdir)/loop_mpi.op-2 -
	mpirun -np 3 loop_mpi | diff $(srcdir)/loop_mpi.op-3 -
	mpirun -np 4 loop_mpi | diff $(srcdir)/loop_mpi.op-4 -

archive:
	rm -rf loop-$(PACKAGE_VERSION)
	mkdir loop-$(PACKAGE_VERSION)
	cp -p $(top_srcdir)/LICENSE.txt $(top_srcdir)/CITATIONS.txt $(srcdir)/README.en $(srcdir)/loop.C $(srcdir)/loop.op $(srcdir)/loop_mpi.C $(srcdir)/loop_mpi.op* $(srcdir)/*.h loop-$(PACKAGE_VERSION)
	tar zcf loop-$(PACKAGE_VERSION).tar.gz loop-$(PACKAGE_VERSION)
	rm -rf loop-$(PACKAGE_VERSION)
