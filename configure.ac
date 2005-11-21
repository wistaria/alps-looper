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
dnl check for PVM library
dnl

AC_MSG_CHECKING([for PVM library])
AC_MSG_RESULT([$ac_cv_alps_have_pvm])
AM_CONDITIONAL(HAVE_PVM, test "$ac_cv_alps_have_pvm" = yes)

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

AC_SUBST([CPPFLAGS_PVM], "$ac_cv_alps_cppflags_pvm")
AC_SUBST([LDFLAGS_PVM], "$ac_cv_alps_ldflags_pvm")
AC_SUBST([LIBS_PVM], "$ac_cv_alps_libs_pvm")

dnl
dnl generate files
dnl

AC_CONFIG_HEADERS(loop_config.h)
AC_CONFIG_FILES(Makefile doc/Makefile test/Makefile ed/Makefile check/Makefile)
AC_CONFIG_FILES(test/model_parameter.ip test/sign.ip ed/diag.ip)
AC_CONFIG_FILES(check/chain-p check/chain-e check/chain.gp check/triangle-p check/triangle-e check/triangle.gp)
AC_OUTPUT
