dnl Process this file with autoconf to produce a configure script.
dnl $Id: configure.ac 705 2004-03-18 01:16:30Z wistaria $

dnl 
dnl Initialization
dnl

dnl disable saving cache
define([AC_CACHE_SAVE], )dnl

AC_PREREQ(2.59)
AC_INIT([looper], 3.0b4)              dnl check here for release!
looper_date="2003/11/20"              dnl check here for release!
AM_INIT_AUTOMAKE

AC_PREFIX_DEFAULT([$HOME/ALPS])

AC_DEFINE_UNQUOTED(LOOPER_VERSION, "$PACKAGE_VERSION",
                   [Define the version of ALPS/looper])
AC_DEFINE_UNQUOTED(LOOPER_DATE, "$looper_date",
                   [Define the published date of ALPS/looper])

AC_MSG_NOTICE([starting configuration for ALPS/looper version $PACKAGE_VERSION ($looper_date)])

dnl
dnl check for ALPS libraries
dnl

AC_MSG_CHECKING([for ALPS Libraries])
if test -f "$ac_cv_alps_include_mk"; then
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

dnl
dnl check for LAPACK
dnl

sinclude(config/ac_lapack.m4)
AC_LAPACK
AM_CONDITIONAL(HAVE_LAPACK, test "$ac_cv_have_lapack" = yes)

dnl
dnl check for MPI library
dnl

AC_MSG_CHECKING([for MPI library])
if test "$ac_cv_alps_have_mpi" = yes; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(HAVE_MPI, test "$ac_cv_alps_have_mpi" = yes)

dnl
dnl check for PVM library
dnl

AC_MSG_CHECKING([for PVM library])
if test "$ac_cv_alps_have_pvm" = yes; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(HAVE_PVM, test "$ac_cv_alps_have_pvm" = yes)

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

dnl
dnl generate files
dnl

AC_CONFIG_HEADERS(looper/config.h)
AC_OUTPUT(Makefile test/Makefile doc/Makefile lattices.xml models.xml percolation-params qmc-params)
