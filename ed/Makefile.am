AM_CPPFLAGS = @CPPFLAGS_SGL@ @LAPACK_CPPFLAGS@ @BINDINGS_CPPFLAGS@
AM_LDFLAGS = -L$(top_builddir) @LDFLAGS_SGL@ @LAPACK_LDFLAGS@
LDADD = @LAPACK_LIBS@ @LIBS_SGL@

check_PROGRAMS = diag

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

diag_SOURCES = diag.C

dist_noinst_SCRIPT = extract.sh run-diag.sh
dist_noinst_DATA = diag.op heisenberg.dat heisenberg.ip ising-x.ip ising-z.ip ising.dat lattices.xml models.xml xy-x.ip xy-z.ip xy.dat

clean: clean-am clean-extra
clean-extra:
	rm -f $(check_PROGRAMS)
distclean: distclean-extra distclean-am
distclean-extra: clean-extra
	rm -rf ./$(DEPDIR)
	rm -f Makefile
