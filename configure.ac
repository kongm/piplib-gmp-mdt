dnl /**-------------------------------------------------------------------**
dnl  **                              PIPLIB                               **
dnl  **-------------------------------------------------------------------**
dnl  **                           configure.ac                            **
dnl  **-------------------------------------------------------------------**
dnl  **                   First version: august 11th  2001                **
dnl  **-------------------------------------------------------------------**/
dnl
#############################################################################
## Require for autoconf >= 2.54.
AC_PREREQ([2.54])

#############################################################################
## Fill here the @bug email adress.
AC_INIT([piplib-gmp], [1.3.6-mp], [cedric.bastoul@inria.fr])
## A common file, which should be included by every source file.
AC_CONFIG_SRCDIR([include/piplib-gmp/piplib-gmp.h])
## Put as most as possible configuration files to an auxialiry
## directory.
AC_CONFIG_AUX_DIR([config])
## Initialize automake. Here, a special tar version that enables
## (very) long filenames.
AM_INIT_AUTOMAKE([1.9 tar-ustar no-define foreign dist-bzip2])
#############################################################################

## Use C language.
AC_PROG_CC
## Use libtool to generate libraries.
AC_PROG_LIBTOOL


## Check for some essential headers.
AC_HEADER_STDC
AC_CHECK_HEADERS([errno.h stddef.h stdlib.h string.h strings.h unistd.h])

## Check for some essential C properties.
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_C_VOLATILE

## Selection of the list of libraries created (bit version).

## Conditional compilation of binaries and libraries.
AC_ARG_WITH(pip,
	AC_HELP_STRING([--with-pip], [Compile PIP binaries]),
        [with_pip=$withval], [with_pip=yes])
AM_CONDITIONAL([WANT_PIP], [test x$with_pip = xyes])

AC_ARG_WITH(libs,
	AC_HELP_STRING([--with-libs], [Compile PIP libraries]),
        [with_libs=$withval], [with_libs=yes])
AM_CONDITIONAL([WANT_LIBS], [test x$with_libs = xyes])

## GMP options.
AC_ARG_WITH(gmp,
        [  --with-gmp=DIR          DIR where the gmp package is installed],
	[ echo "Package gmp : $withval";
	  CPPFLAGS="${CPPFLAGS} -I$withval/include";
	  LDFLAGS="${LDFLAGS} -L$withval/lib"
          LIBS="${LIBS} -lgmp"
	  AM_CONDITIONAL([WANT_LIBPIPLIBMP], [true])])

AC_ARG_WITH(gmp_include,
        [  --with-gmp-include=DIR  DIR where gmp.h is installed],
	[ echo "Package gmp-include : $withval";
	  CPPFLAGS="${CPPFLAGS} -I$withval";
	  AM_CONDITIONAL([WANT_LIBPIPLIBMP], [true])])

AC_ARG_WITH(gmp_library,
        [  --with-gmp-library=DIR  DIR where the gmp library is installed],
	[ echo "Package gmp-library : $withval";
	  LDFLAGS="${LDFLAGS} -L$withval"
          LIBS="${LIBS} -lgmp"
	  AM_CONDITIONAL([WANT_LIBPIPLIBMP], [true])])

#############################################################################
## Variable substitutions.
AC_SUBST(ac_aux_dir)

#############################################################################
## Configure Makefiles.
AC_CONFIG_FILES([
	Makefile
	doc/Makefile
	example/Makefile
	include/Makefile
	source/Makefile
	test/Makefile
	],
	[test -z "$CONFIG_HEADERS" || echo timestamp > source/stamp-h.in])

## Output the result.
AC_OUTPUT

echo "             /*-----------------------------------------------*"
echo "              *         PIP/PipLib-GMP configuration is OK    *"
echo "              *-----------------------------------------------*/"
echo "It appears that your system is OK to start PIP/PipLib compilation. You"
echo "need now to type \"make uninstall\" if your are upgrading an old version."
echo "Next, type \"make\" to build the program. Lastly type \"make install\" to"
echo "install it on your system (log as root if necessary)."
