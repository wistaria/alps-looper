dnl Process this file with autoconf to produce a configure script.

dnl disable saving cache
define([AC_CACHE_SAVE], )dnl
AC_PREREQ(2.59)

dnl
dnl Version Section ---------------------------------------------------------
dnl

dnl name of package
define([package_name], [looper])

dnl version of package
define([package_version], [esyscmd(config/extract_version.sh)])
define([package_date], [esyscmd(config/extract_date.sh)])

dnl 
dnl Initialization ----------------------------------------------------------
dnl

AC_INIT(package_name, package_version)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT([$HOME/ALPS])
AM_MAINTAINER_MODE

AC_MSG_NOTICE([starting configuration for ALPS/looper version package_version (package_date)])

undefine([package_name])
undefine([package_version])
undefine([package_date])

dnl
dnl check for ALPS library --------------------------------------------------
dnl

AC_MSG_CHECKING([for ALPS Libraries])
if test -f "$ac_cv_alps_includedir/alps/config.h"; then
  AC_MSG_RESULT([OK])
else
  AC_MSG_RESULT([not found])
  AC_MSG_ERROR([ALPS Libraries not found.])
fi

AC_MSG_CHECKING([for ALPS version])
AC_MSG_RESULT([$ac_cv_alps_version])

if test "$ac_cv_alps_type" = "full"; then :; else
  AC_MSG_NOTICE([detected ALPS Light Libraries])
  AC_MSG_ERROR([The full version of ALPS Libraries is required.])
fi

if test -z "${enable_shared+set}" && test -n "$ac_cv_enable_shared"; then
  enable_shared="$ac_cv_enable_shared"
fi
if test -z "${enable_static+set}" && test -n "$ac_cv_enable_static"; then
  enable_static="$ac_cv_enable_static"
fi

AC_PROG_CC
AC_PROG_CXX
F77=no
AC_PROG_LIBTOOL
AC_LANG_CPLUSPLUS

CPPFLAGS="$ac_cv_alps_cppflags"
CXXFLAGS="$ac_cv_alps_cxxflags"
LDFLAGS="$ac_cv_alps_ldflags"
LIBS="$ac_cv_alps_libs"

dnl
dnl check for MPI library
dnl

AC_MSG_CHECKING([for MPI library])
AC_MSG_RESULT([$ac_cv_alps_have_mpi])
AM_CONDITIONAL(HAVE_MPI, test "$ac_cv_alps_have_mpi" = yes)

dnl
dnl check for LAPACK
dnl

AC_MSG_CHECKING([for LAPACK library])
AC_MSG_RESULT([$ac_cv_have_lapack])
if test "$ac_cv_have_lapack" != yes; then
  AC_MSG_ERROR([LAPACK Library not found.])
fi

dnl
dnl compile flags
dnl

CXXFLAGS="$ac_cv_alps_cxxflags -DBOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS"
CPPFLAGS="$ac_cv_alps_cppflags"
LDFLAGS="$ac_cv_alps_ldflags"
LIBS=
AC_SUBST([DEFAULT_INCLUDES], '-I$(top_srcdir) -I$(top_builddir)')
AC_SUBST([LDADD], "$ac_cv_alps_libs")

AC_SUBST([CPPFLAGS_MPI], "$ac_cv_alps_cppflags_mpi")
AC_SUBST([LDFLAGS_MPI], "$ac_cv_alps_ldflags_mpi")
AC_SUBST([LIBS_MPI], "$ac_cv_alps_libs_mpi")

dnl
dnl methods
dnl

AC_ARG_ENABLE(all-method,
  AC_HELP_STRING([--enable-all-method],
    [enable all methods]),
  [
  if test "x$enableval" != "xno"; then
    enable_diag=yes
    enable_ising=yes
  else
    enable_pi=no
    enable_sse=no
    enable_sse_qwl=no
  fi
  ]
)

AC_ARG_ENABLE(diag,
  AC_HELP_STRING([--enable-diag],
    [enable exact diagonalization method]),
  [test "x$enableval" != "xno" && enable_diag=yes]
)
test -z "$enable_diag" && enable_diag=no
AM_CONDITIONAL(ENABLE_DIAG, test "$enable_diag" = yes)

AC_ARG_ENABLE(ising,
  AC_HELP_STRING([--enable-ising],
    [enable Swendsen-Wang method for classical Ising model]),
  [test "x$enableval" != "xno" && enable_ising=yes]
)
test -z "$enable_ising" && enable_ising=no
AM_CONDITIONAL(ENABLE_ISING, test "$enable_ising" = yes)

AC_ARG_ENABLE(path-integral,
  AC_HELP_STRING([--disable-path-integral],
    [disable path integral method]),
  [test "x$enableval" != "xyes" && enable_pi=no]
)
test -z "$enable_pi" && enable_pi=yes
AM_CONDITIONAL(ENABLE_PI, test "$enable_pi" = yes)

AC_ARG_ENABLE(sse,
  AC_HELP_STRING([--disable-sse],
    [disable SSE method]),
  [test "x$enableval" != "xyes" && enable_sse=no]
)
test -z "$enable_sse" && enable_sse=yes
AM_CONDITIONAL(ENABLE_SSE, test "$enable_sse" = yes)

AC_ARG_ENABLE(sse-qwl,
  AC_HELP_STRING([--disable-sse-qwl],
    [disable SSE quantum Wang-Landau method]),
  [test "x$enableval" != "xyes" && enable_sse_qwl=no]
)
test -z "$enable_sse_qwl" && enable_sse_qwl=yes
AM_CONDITIONAL(ENABLE_SSE_QWL, test "$enable_sse_qwl" = yes)

dnl
dnl subdirectories
dnl

if test -f "$srcdir/check/Makefile.in"; then
  exist_check=yes
fi
AM_CONDITIONAL(EXIST_CHECK, test "$exist_check" = yes)
AM_CONDITIONAL(BUILD_CHECK, false)

if test -f "$srcdir/extras/Makefile.in"; then
  exist_extras=yes
fi
AM_CONDITIONAL(EXIST_EXTRAS, test "$exist_extras" = yes)
AM_CONDITIONAL(BUILD_EXTRAS, false)

if test -f "$srcdir/test/Makefile.in"; then
  exist_test=yes
fi
AM_CONDITIONAL(EXIST_TEST, test "$exist_test" = yes)
AM_CONDITIONAL(BUILD_TEST, false)

if test -f "$srcdir/v3.1/Makefile.in"; then
  exist_v31=yes
fi
AM_CONDITIONAL(EXIST_V31, test "$exist_v31" = yes)
AM_CONDITIONAL(BUILD_V31, false)

dnl
dnl generate files
dnl

AC_CONFIG_FILES(Makefile doc/Makefile)

if test "$exist_check" = yes; then
  AC_CONFIG_FILES(check/Makefile check/site.gp check/chain.gp check/triangle.gp)
fi
if test "$exist_extras" = yes; then
  AC_CONFIG_FILES(extras/Makefile extras/gap/Makefile)
fi
if test "$exist_test" = yes; then
  AC_CONFIG_FILES(test/Makefile)
fi
if test "$exist_v31" = yes; then
  AC_CONFIG_FILES(v3.1/Makefile)
fi

AC_OUTPUT
