#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(tepla, 2.0, tepla-admin@cipher.risk.tsukuba.ac.jp)

m4_include([m4/ax_check_openssl.m4])

AC_CONFIG_SRCDIR([src/tepla/hash.c])
AC_CONFIG_HEADERS([config.h])

if test "" = "$CFLAGS"; then
   CFLAGS="-Wall -Wextra -pedantic -g -O3 -std=c11"
fi

if test "" = "x$LDFLAGS"; then
   LDFLAGS=""
fi

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET

# Checks for gmp
AC_CHECK_LIB([gmp], [main], [], [AC_MSG_ERROR("Test for library gmp failed.")])
AC_CHECK_HEADER([gmp.h], [], [AC_MSG_ERROR("cannot find gmp.h")])

# Checks for openssl
#AC_CHECK_LIB([crypto], [main], [], [AC_MSG_ERROR("Test for library openssl crypto failed.")])
#AC_CHECK_HEADER([openssl/sha.h], [], [AC_MSG_ERROR("cannot find openssl/sha.h")])
AX_CHECK_OPENSSL([],[AC_MSG_ERROR("Test for library openssl failed")])

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

# Select algorithm in Library
AC_ARG_ENABLE([fastalg], [ --enable-fastalg turn on to use faster algorithms [default yes]],[],[enable_fastalg=yes])
AC_MSG_CHECKING(whether to enable fastalg)
if test "x$enable_fastalg" = "xyes"; then
   AC_MSG_RESULT(yes)
   AC_DEFINE(ENABLE_FASTALG,[], enable faster algorithms)
else
   AC_MSG_RESULT(no)
fi

# Select Library type (static or shared)
AC_ARG_ENABLE([shared], [ --enable-shared turn on to make shared library [default no]],[],[enable_shared=no])
AC_MSG_CHECKING(whether to enable shared)
if test "x$enable_shared" = "xyes"; then
   AC_MSG_RESULT(yes)
   AC_SUBST(LIBTYPE,"shared")
else
   AC_MSG_RESULT(no)
   AC_SUBST(LIBTYPE,"static")
fi

AC_CONFIG_FILES([
		Makefile
		src/tepla/Makefile
		src/tepla/ec_bn254/Makefile
		src/tepla/ec_bn254/test/Makefile])

AC_OUTPUT
