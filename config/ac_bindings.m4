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
        for s in bindings ublas_bindings boost-sandbox bindings_2003_09_11
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
    BINDINGS_CPPFLAGS="-I$bindings_dir -DBOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK"
    AC_DEFINE(HAVE_BINDINGS, [], [Define if you have uBLAS LAPACK bindings.])
  else
    ac_cv_have_bindings=no
    ac_cv_bindings_dir=
    BINDINGS_CPPFLAGS=
  fi
  ]
)
