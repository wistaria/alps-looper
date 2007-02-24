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
define([package_bugreport], [wistaria@comp-phys.org])

dnl 
dnl Initialization ----------------------------------------------------------
dnl

AC_INIT(package_name, package_version, package_bugreport)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT([$HOME/ALPS])
AM_MAINTAINER_MODE

if test -f "$srcdir/looper/gap.h"; then
  looper_pro=yes
  AC_MSG_NOTICE([starting configuration for ALPS/looper professional version package_version (package_date)])
else
  AC_MSG_NOTICE([starting configuration for ALPS/looper version package_version (package_date)])
fi
AM_CONDITIONAL(LOOPER_PRO, test "$looper_pro" = "yes")

undefine([package_name])
undefine([package_version])
undefine([package_date])
undefine([package_bugreport])

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
dnl check for ALPS/parapack
dnl

AC_CHECK_HEADER([parapack/version.h])
test "$ac_cv_header_parapack_version_h" = yes && ac_cv_have_parapack=yes
test "$ac_cv_have_parapack" = yes && AC_DEFINE([HAVE_PARAPACK])
AM_CONDITIONAL(HAVE_PARAPACK, test "$ac_cv_have_parapack" = yes)

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
dnl config file
dnl

config_file=
AC_ARG_WITH(config,
  AC_HELP_STRING([--with-config=FILE],
    [use config file instead of the default one (loop_config.h)]),
  [
  if test "x$withval" != xno; then
    config_file=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
  fi
  ]
)
test -n "$config_file" && AC_MSG_NOTICE([using config file: $config_file])

dnl
dnl compile flags
dnl

CXXFLAGS="$ac_cv_alps_cxxflags -DBOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS"
CPPFLAGS="$ac_cv_alps_cppflags"
LDFLAGS="$ac_cv_alps_ldflags"
LIBS=
AC_SUBST([DEFAULT_INCLUDES], '-I$(top_srcdir) -I$(top_builddir)')
if test -n "$config_file"; then
  DEFAULT_INCLUDES="$DEFAULT_INCLUDES -DLOOP_CONFIG_HEADER=\\\"$config_file\\\""
fi
AC_SUBST([LDADD], "$ac_cv_alps_libs")

AC_SUBST([LDFLAGS_MPI], "$ac_cv_alps_ldflags_mpi")
AC_SUBST([LIBS_MPI], "$ac_cv_alps_libs_mpi")

AC_SUBST([LDFLAGS_PS], "$ac_cv_alps_ldflags_mpi")
AC_SUBST([LIBS_PS], "-lparapack_sgl $ac_cv_alps_libs_mpi")

AC_SUBST([LDFLAGS_PP], "$ac_cv_alps_ldflags_mpi")
AC_SUBST([LIBS_PP], "-lparapack_mpi $ac_cv_alps_libs_mpi")

dnl
dnl methods
dnl

AC_ARG_ENABLE(diag,
  AC_HELP_STRING([--disable-diag],
    [disable exact diagonalization method]),
  [test "x$enableval" != "xyes" && enable_diag=no]
)
test -f "$srcdir/diag.C" || enable_diag=no
test -z "$enable_diag" && enable_diag=yes
AM_CONDITIONAL(ENABLE_DIAG, test "$enable_diag" = yes)

AC_ARG_ENABLE(ising,
  AC_HELP_STRING([--disable-ising],
    [disable Swendsen-Wang method for classical Ising model]),
  [test "x$enableval" != "xyes" && enable_ising=no]
)
test -f "$srcdir/ising.C" || enable_ising=no
test -z "$enable_ising" && enable_ising=yes
AM_CONDITIONAL(ENABLE_ISING, test "$enable_ising" = yes)

AC_ARG_ENABLE(path-integral,
  AC_HELP_STRING([--disable-path-integral],
    [disable path integral method]),
  [test "x$enableval" != "xyes" && enable_pi=no]
)
test -f "$srcdir/path_integral.C" || enable_pi=no
test -z "$enable_pi" && enable_pi=yes
AM_CONDITIONAL(ENABLE_PI, test "$enable_pi" = yes)

AC_ARG_ENABLE(sse,
  AC_HELP_STRING([--disable-sse],
    [disable SSE method]),
  [test "x$enableval" != "xyes" && enable_sse=no]
)
test -f "$srcdir/sse.C" || enable_sse=no
test -z "$enable_sse" && enable_sse=yes
AM_CONDITIONAL(ENABLE_SSE, test "$enable_sse" = yes)

AC_ARG_ENABLE(sse-qwl,
  AC_HELP_STRING([--disable-sse-qwl],
    [disable SSE quantum Wang-Landau method]),
  [test "x$enableval" != "xyes" && enable_sse_qwl=no]
)
test -f "$srcdir/sse_qwl.C" || enable_sse_qwl=no
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
  AC_CONFIG_FILES(extras/Makefile)
  AC_CONFIG_FILES(extras/gap/Makefile)
  AC_CONFIG_FILES(extras/localsus/Makefile)
  AC_CONFIG_FILES(extras/top/Makefile)
  AC_CONFIG_FILES(extras/transmag/Makefile extras/transmag/alternating_chain.xml)
fi
if test "$exist_test" = yes; then
  AC_CONFIG_FILES(test/Makefile)
fi
if test "$exist_v31" = yes; then
  AC_CONFIG_FILES(v3.1/Makefile)
fi

AC_OUTPUT
