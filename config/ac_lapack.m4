AC_DEFUN([AC_LAPACK],
  [
  AC_SUBST(LAPACK_CPPFLAGS)
  AC_SUBST(LAPACK_LDFLAGS)
  AC_SUBST(LAPACK_LIBS)
  
  AC_ARG_WITH(lapack,
    AC_HELP_STRING([--with-lapack=DIR],[use LAPACK Library]),
    [
    if test "x$withval" = "xno"; then
      lapack=no
    else
      lapack=yes
      if test "x$withval" != "xyes"; then
        lapack_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  AC_ARG_WITH(mkl,
    AC_HELP_STRING([--with-mkl=PROC],
    [use Intel Math Kernel Library (PROC=p3, p4, itp)]),
    [
    case "x$withval" in
      xno)
        mkl=no
        ;;
      xp3)
        mkl=yes
	mkl_proc="p3"
        ;;
      xp4)
        mkl=yes
	mkl_proc="p4"
        ;;
      xitp)
        mkl=yes
	mkl_proc="itp"
        ;;
      *)
        mkl=yes
	mkl_proc="p4"
    esac
    ]
  )

  AC_ARG_WITH(scsl,
    AC_HELP_STRING([--with-scsl],
    [use SGI SCSL Scientific Library]),
    [
    case "x$withval" in
      xno)
        scsl=no
        ;;
      *)
        scsl=yes
    esac
    ]
  )

  if test "$mkl" != no; then

    AC_MSG_CHECKING([for Math Kernel Library])
    if test -n "$lapack_dir" && test -d "$lapack_dir"; then
      mkl_basedir="$lapack_dir"
      AC_MSG_RESULT([$lapack_dir])
    elif test -d "/opt/intel/mkl61"; then
      mkl_basedir=/opt/intel/mkl61
      AC_MSG_RESULT([$mkl_basedir])
    elif test -d "/opt/intel/mkl"; then
      mkl_basedir=/opt/intel/mkl
      AC_MSG_RESULT([$mkl_basedir])
    else
      AC_MSG_RESULT([OK])
    fi

    AC_MSG_CHECKING([for processor type])
    case "x$mkl_proc" in
      xp3)
        mkl_procname="Pentium 3"
        mkl_dir="32"
        mkl_lib="mkl_p3"
        ;;
      xp4)
        mkl_procname="Pentium 4"
        mkl_dir="32"
        mkl_lib="mkl_p4"
        ;;
      xitp)
        mkl_procname="Itanium"
        mkl_dir="64"
        mkl_lib="mkl_itp"
        ;;
      *)
        mkl_procname="Pentium 4"
        mkl_dir="32"
        mkl_lib="mkl_p4"
    esac
    if test "$mkl_dir" = "32"; then
      if test -f "$mkl_basedir/lib/$mkl_dir/libmkl_ia32.a"; then
        # version 6
        mkl_lib="mkl_ia32"
      fi
    fi
    if test "$mkl_dir" = "64"; then
      if test -f "$mkl_basedir/lib/$mkl_dir/libmkl_ipf.a"; then
        # version 6
        mkl_lib="mkl_ipf"
      fi
    fi
    AC_MSG_RESULT([$mkl_procname])
    LAPACK_LDFLAGS="-L$mkl_basedir/lib/$mkl_dir"

    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$LAPACK_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$LAPACK_LDFLAGS $LDFLAGS"
  
    found=no
  
    LAPACK_LIBS=""
    if test "$mkl_proc" != "itp"; then
      AC_CHECK_FUNC(d_abs,,
      [
      AC_CHECK_LIB(g2c, d_abs, 
      [
      LAPACK_LIBS="-lg2c $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ],
      [
      AC_CHECK_LIB(F90, d_abs, 
      [
      LAPACK_LIBS="-lF90 $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ])
      ])
      ])
      AC_CHECK_LIB(guide, omp_in_parallel_, 
      [
      LAPACK_LIBS="-lguide $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      AC_CHECK_LIB($mkl_lib, sasum_, 
      [
      LAPACK_LIBS="-l$mkl_lib $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      AC_CHECK_LIB(mkl_lapack, dsyev_,
        [LAPACK_LIBS="-lmkl_lapack $LAPACK_LIBS"; found=yes])
      ])
      ])
    else
      AC_CHECK_FUNC(d_abs,,
      [
      AC_CHECK_LIB(g2c, d_abs, 
      [
      LAPACK_LIBS="-lg2c $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ],
      [
      AC_CHECK_LIB(F90, d_abs, 
      [
      LAPACK_LIBS="-lF90 $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ])
      ])
      ])
      AC_CHECK_FUNC(f_concat,,
      [
      AC_CHECK_LIB(IEPCF90, f_concat, 
      [
      LAPACK_LIBS="-lIEPCF90 $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ])
      ])
      AC_CHECK_LIB(guide, omp_in_parallel_, 
      [
      LAPACK_LIBS="-lguide $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      AC_CHECK_LIB($mkl_lib, sasum_, 
      [
      LAPACK_LIBS="-l$mkl_lib $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      AC_CHECK_LIB(mkl_lapack, dsyev_,
        [LAPACK_LIBS="-lmkl_lapack $LAPACK_LIBS"; found=yes])
      ])
      ])
    fi

    if test "$found" = no; then
      if test "$mkl" = yes; then
        AC_MSG_ERROR([check for Math Kernel Library failed])
      else
        mkl=no
      fi
    else
      mkl=yes
    fi

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE
  fi

  if test "$mkl" != yes && test "$scsl" != no; then
    AC_MSG_NOTICE([checking for SGI SCSL library])
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$LAPACK_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$LAPACK_LDFLAGS $LDFLAGS"

    found=no

    if test "$found" = no; then
      AC_CHECK_LIB(scs, dsyev_, [LAPACK_LIBS="-lscs"; found=yes])
    fi

    if test "$found" = no; then
      if test "$scsl" = yes; then
        AC_MSG_ERROR([check for SGI SCSL library failed])
      else
        scsl=no
      fi
    else
      scsl=yes
    fi

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE
  fi

  if test "$mkl" != yes && test "$scsl" != yes && test "$lapack" != "no"; then
    AC_MSG_CHECKING([for LAPACK root directory])
    if test "x$lapack_dir" = "x"; then
      for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src; do
        if test -f "$d/lib/liblapack.a"; then
          lapack_dir="$d"
          break
        fi
        if test -f "$d/lapack/lib"; then
          lapack_dir="$d/lapack"
          break
        fi
      done
      if test -n "$lapack_dir"; then
        AC_MSG_RESULT([$lapack_dir])
      else
        AC_MSG_RESULT([OK])
      fi
    else
      AC_MSG_RESULT([$lapack_dir])
      if test -d "$lapack_dir"; then :; else
        AC_MSG_ERROR([$lapack_dir not found])
      fi
    fi

    if test -n "$lapack_dir"; then
      if test -f "$lapack_dir/lib"; then
        LAPACK_LDFLAGS="-L$lapack_dir/lib"
      fi
    fi

    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$LAPACK_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$LAPACK_LDFLAGS $LDFLAGS"

    found=no

    if test "$found" = no; then
      AC_CHECK_FUNC(d_abs,,
      [
      AC_CHECK_LIB(g2c, d_abs, 
      [
      LAPACK_LIBS="-lg2c $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ],
      [
      AC_CHECK_LIB(F90, d_abs, 
      [
      LAPACK_LIBS="-lF90 $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      ])
      ])
      ])
      AC_CHECK_LIB(blas, sasum_, 
      [
      LAPACK_LIBS="-lblas $LAPACK_LIBS"
      LIBS="$LAPACK_LIBS $ac_save_LIBS"
      AC_CHECK_LIB(lapack, dsyev_,
        [LAPACK_LIBS="-llapack $LAPACK_LIBS"; found=yes])
      ])
    fi

    # for Digital Extended Math Library
    if test "$found" = no; then
      AC_MSG_NOTICE([checking for Digital Extended Math Library])
      LAPACK_CPPFLAGS=
      LAPACK_LDFLAGS=
      LAPACK_LIBS=
      CPPFLAGS="$ac_save_CPPFLAGS"
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_CHECK_LIB(dxml, sasum_, 
      [
      AC_CHECK_LIB(dxml, dsyev_, [LAPACK_LIBS="-ldxml"; found=yes])
      ])
    fi

    # for vecLib on Mac OS X
    if test "$found" = no; then
      AC_MSG_NOTICE([checking for vecLib framework on Mac OS X])
      LAPACK_CPPFLAGS=
      LAPACK_LDFLAGS=
      LAPACK_LIBS=
      CPPFLAGS="$ac_save_CPPFLAGS -faltivec -framework vecLib"
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_CHECK_FUNC(dsyev_,
      [
      LAPACK_CPPFLAGS="-faltivec -framework vecLib"; found=yes
      ])
    fi

    if test "$found" = no; then
      if test "$lapack" = yes; then
        AC_MSG_ERROR([check for LAPACK library failed])
      fi
      lapack=no
    fi

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE
  fi
  
  if test "$lapack" != no; then
    AC_MSG_NOTICE([LAPACK enabled])
    AC_DEFINE(HAVE_LAPACK)
    ac_cv_have_lapack=yes
  else
    AC_MSG_NOTICE([LAPACK disabled])
    ac_cv_have_lapack=no
    LAPACK_CPPFLAGS=
    LAPACK_LDFLAGS=
    LAPACK_LIBS=
  fi
  ]
)
