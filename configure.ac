dnl configuration script for ALPS/looper
dnl Process this file with autoconf to produce a configure script.

dnl disable saving cache
define([AC_CACHE_SAVE], )dnl
AC_PREREQ(2.59)

dnl
dnl Version Section ---------------------------------------------------------
dnl

define([looper_version], [esyscmd(config/extract_version.sh)])
define([looper_date], [esyscmd(config/extract_date.sh)])

dnl 
dnl Initialization
dnl

AC_INIT([looper], lopper_version)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT([$HOME/ALPS])

AC_MSG_NOTICE([starting configuration for ALPS/looper version looper_version (looper_date)])

undefine([looper_version])
undefine([looper_date])

dnl
dnl check for ALPS libraries
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
  AC_MSG_ERROR([ALPS/looper requies the full version of ALPS Libraries])
fi

AC_PROG_CXX
AC_PROG_LIBTOOL

CPPFLAGS="$ac_cv_alps_cppflags"
CXXFLAGS="$ac_cv_alps_cxxflags"
LDFLAGS="$ac_cv_alps_ldflags"
LIBS="$ac_cv_alps_libs"

AC_SUBST([CXXFLAGS], $ac_cv_alps_cxxflags)
AC_SUBST([DEFAULT_INCLUDES], '-I$(top_srcdir) -I$(top_builddir)')

AC_SUBST([CPPFLAGS_SGL], $ac_cv_alps_cppflags)
AC_SUBST([LDFLAGS_SGL], $ac_cv_alps_ldflags)
AC_SUBST([LIBS_SGL], $ac_cv_alps_libs)

AC_SUBST([CPPFLAGS_MPI], $ac_cv_alps_cppflags_mpi)
AC_SUBST([LDFLAGS_MPI], $ac_cv_alps_ldflags_mpi)
AC_SUBST([LIBS_MPI], $ac_cv_alps_libs_mpi)

AC_SUBST([CPPFLAGS_PVM], $ac_cv_alps_cppflags_pvm)
AC_SUBST([LDFLAGS_PVM], $ac_cv_alps_ldflags_pvm)
AC_SUBST([LIBS_PVM], $ac_cv_alps_libs_pvm)

AC_LANG_CPLUSPLUS

dnl
dnl check for LAPACK
dnl

AC_LAPACK
AM_CONDITIONAL(HAVE_LAPACK, test "$ac_cv_alps_have_lapack" = yes)
if test "$ac_cv_have_lapack" != yes; then
  AC_MSG_ERROR([LAPACK Library not found.])
fi

dnl
dnl check for uBLAS LAPACK bindings
dnl

AC_BINDINGS
if test "$ac_cv_have_bindings" != yes; then
  AC_MSG_ERROR([uBLAS LAPACK bindings not found.])
fi

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
dnl generate files
dnl

AC_CONFIG_HEADERS(looper/config.h)
AC_CONFIG_FILES(Makefile doc/Makefile test/Makefile ed/Makefile)
AC_CONFIG_FILES(test/xxz_model.ip ed/diag.ip)
AC_OUTPUT
