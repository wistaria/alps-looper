AM_CPPFLAGS = @CPPFLAGS_SGL@ @LAPACK_CPPFLAGS@
AM_LDFLAGS = -L$(top_builddir) @LDFLAGS_SGL@ @LAPACK_LDFLAGS@
LDADD = @LAPACK_LIBS@ @LIBS_SGL@

check_PROGRAMS =
if HAVE_LAPACK
check_PROGRAMS += diag
endif

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/run-test $(srcdir)

diag_SOURCES = diag.C
