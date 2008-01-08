dist_noinst_DATA = loop.op
check_PROGRAMS = loop
loop_SOURCES = loop.C common.h observable.h options.h union_find.h

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
