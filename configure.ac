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
dnl check for ALPS and ALPS/parapack
dnl

AC_ALPS
AC_PARAPACK

dnl
dnl compiler
dnl

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

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

AC_SUBST([DEFAULT_INCLUDES], '-I$(top_srcdir) -I$(top_builddir)')
if test -n "$config_file"; then
  DEFAULT_INCLUDES="$DEFAULT_INCLUDES -DLOOP_CONFIG_HEADER=\\\"$config_file\\\""
fi

dnl
dnl algorithms
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
  AC_HELP_STRING([--enable-sse-qwl],
    [enable SSE quantum Wang-Landau method]),
  [test "x$enableval" != "xno" && enable_sse_qwl=yes]
)
test -f "$srcdir/sse_qwl.C" || enable_sse_qwl=no
test -z "$enable_sse_qwl" && enable_sse_qwl=no
AM_CONDITIONAL(ENABLE_SSE_QWL, test "$enable_sse_qwl" = yes)

dnl
dnl subdirectories
dnl

if test -f "$srcdir/check/Makefile.in"; then
  exist_check=yes
fi
AM_CONDITIONAL(EXIST_CHECK, test "$exist_check" = yes)

if test -f "$srcdir/extras/Makefile.in"; then
  exist_extras=yes
fi
AM_CONDITIONAL(EXIST_EXTRAS, test "$exist_extras" = yes)

if test -f "$srcdir/simple/Makefile.in"; then
  exist_simple=yes
fi
AM_CONDITIONAL(EXIST_SIMPLE, test "$exist_simple" = yes)

if test -f "$srcdir/test/Makefile.in"; then
  exist_test=yes
fi
AM_CONDITIONAL(EXIST_TEST, test "$exist_test" = yes)

if test -f "$srcdir/v3.1/Makefile.in"; then
  exist_v31=yes
fi
AM_CONDITIONAL(EXIST_V31, test "$exist_v31" = yes)
AM_CONDITIONAL(BUILD_V31, false)

dnl
dnl compiler options
dnl

CPPFLAGS="$CPPFLAGS -DBOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS"

dnl
dnl generate files
dnl

AC_CONFIG_FILES(Makefile doc/Makefile)

if test "$exist_check" = yes; then
  AC_CONFIG_FILES(check/Makefile check/chain.gp check/dterm.gp check/kagome.gp check/plot.gp check/site.gp check/transmag.gp check/triangle.gp check/xyz.gp)
fi
if test "$exist_extras" = yes; then
  AC_CONFIG_FILES(extras/Makefile)
  AC_CONFIG_FILES(extras/corrlen/Makefile)
  AC_CONFIG_FILES(extras/gap/Makefile)
  AC_CONFIG_FILES(extras/localsus/Makefile)
  AC_CONFIG_FILES(extras/qwl/Makefile)
  AC_CONFIG_FILES(extras/top/Makefile)
  AC_CONFIG_FILES(extras/transmag/Makefile extras/transmag/alternating_chain.xml)
  AC_CONFIG_FILES(extras/xyz/Makefile extras/xyz/xyz_model.xml)
fi
if test "$exist_simple" = yes; then
  AC_CONFIG_FILES(simple/Makefile)
fi
if test "$exist_test" = yes; then
  AC_CONFIG_FILES(test/Makefile)
fi
if test "$exist_v31" = yes; then
  AC_CONFIG_FILES(v3.1/Makefile)
fi

AC_OUTPUT
