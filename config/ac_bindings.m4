#
# check for Boost uBLAS lapack bindings
#

AC_DEFUN([AC_BINDINGS],
  [
  AC_SUBST(BINDINGS_CPPFLAGS)
    
  AC_ARG_WITH(bindings,
    AC_HELP_STRING([--with-bindings=DIR],
                   [path to Boost uBLAS bindings root directory]),
    [
    if test "x$withval" = xno; then
      bindings=no
    else
      bindings=yes
      if test "x$withval" != xyes; then
        bindings_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  found=no

  if test "x$bindings" != xno; then
    AC_MSG_CHECKING([for Boost uBLAS bindings library])
    if test -n "$bindings_dir"; then
      if test -f "$bindings_dir/boost/numeric/bindings/lapack/lapack.hpp"; then
        found=yes
      else
        bindings_dir=no
      fi
    else
      bindings_dir=no
      for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src
      do
        for s in boost-sandbox boost-bindings bindings ublas_bindings
        do
          if test -f "$d/$s/boost/numeric/bindings/lapack/lapack.hpp"; then
            found=yes
            bindings_dir="$d/$s"
            break
          fi
        done
        if test "$found" = yes; then
          break
        fi
      done
    fi
    AC_MSG_RESULT([$bindings_dir])
  fi

  if test "$found" = yes; then
    ac_cv_have_bindings=yes
    ac_cv_bindings_dir="$bindings_dir"
    BINDINGS_CPPFLAGS="-I$bindings_dir -DBOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS"
    AC_DEFINE(HAVE_BINDINGS, [], [Define if you have uBLAS LAPACK bindings.])
  else
    ac_cv_have_bindings=no
    ac_cv_bindings_dir=
    BINDINGS_CPPFLAGS=
  fi

  if test "$ac_cv_have_bindings" = yes; then
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $BINDINGS_CPPFLAGS"
    AC_CHECK_HEADER([boost/numeric/bindings/lapack/syev.hpp])
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_LANG_RESTORE
  
    if test "$ac_cv_header_boost_numeric_bindings_lapack_syev_hpp" = yes; then
      AC_DEFINE(HAVE_BOOST_NUMERIC_BINDINGS_LAPACK_SYEV_HPP)
    fi
  fi
  ]
)
