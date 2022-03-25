AC_INIT([hgwrr], 0.1-0)

AC_LANG(C++)
AC_REQUIRE_CPP

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test "${GSL_CONFIG}" != ""; then
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs-without-cblas`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_OUTPUT(src/Makevars)