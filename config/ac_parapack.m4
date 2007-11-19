dnl check for ALPS/parapack
AC_DEFUN([AC_PARAPACK],
  [
  AC_CHECK_HEADER([parapack/version.h])
  test "$ac_cv_header_parapack_version_h" = yes && ac_cv_have_parapack=yes
  test "$ac_cv_have_parapack" = yes && AC_DEFINE([HAVE_PARAPACK])
  AM_CONDITIONAL(HAVE_PARAPACK, test "$ac_cv_have_parapack" = yes)
  if test "$ac_cv_have_parapack" = yes; then
    LDADD="-lparapack_sgl $LDADD"
    if test "$ac_cv_have_mpi" = yes; then
      LIBS_MPI="-lparapack_mpi $LIBS_MPI"
    fi
  fi
  ]
)
