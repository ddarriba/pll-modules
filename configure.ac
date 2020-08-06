#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([pll-modules], [0.1.0], [Diego.Darriba@h-its.org])
AM_INIT_AUTOMAKE([subdir-objects])
AC_LANG([C])
AC_CONFIG_SRCDIR([src/pllmod_common.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

# set default gcc options
CFLAGS='-g'

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_INSTALL

if test "x$LEX" != xflex; then
  AC_MSG_ERROR(could not find required installation of FLEX)
fi

if test "x$YACC" != x"bison -y"; then
  AC_MSG_ERROR(could not find required installation of BISON)
fi

LT_INIT


# Checks for libraries.
#PLL_VERSION=0.2.0
#AM_PATH_PLL($PLL_VERSION,
#            :,
#            AC_MSG_ERROR([*** PLL version $PLL_VERSION not found!])
#)
#EXTRA_CFLAGS="$EXTRA_CFLAGS $PLL_CFLAGS"
#EXTRA_LDFLAGS="$EXTRA_LDFLAGS $PLL_LIBS"

AC_CHECK_LIB([m],[exp])

# Checks for header files.
AC_CHECK_HEADERS([assert.h math.h stdio.h stdlib.h string.h ctype.h x86intrin.h])
AC_CHECK_HEADERS([pll.h], [], [AC_MSG_ERROR([pll.h not found])])
#PKG_CHECK_MODULES([PLL], [libpll], [have_pll=yes], [have_pll=no])
AM_CONDITIONAL(HAVE_PLL_DPKG, test "x${have_pll}" = "xyes")

#AC_CHECK_HEADERS([pll.h],
#        [mypj_found_int_headers=yes; break;])

#AS_IF([test "x$mypj_found_int_headers" != "xyes"],
#        [AC_MSG_ERROR([Unable to find PLL header])])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([malloc realloc asprintf memcpy memset])

have_avx=no
have_sse=no
have_ps2pdf=no

AX_EXT

AS_IF([test "x${ax_cv_have_avx_ext}" = "xyes"], [
AC_ARG_ENABLE(avx, AS_HELP_STRING([--enable-avx], [Build with AVX support]), [have_avx=yes])
])

AS_IF([test "x${ax_cv_have_sse3_ext}" = "xyes"], [
AC_ARG_ENABLE(sse, AS_HELP_STRING([--enable-sse],[Build with SSE support]), [have_sse=yes])
])


AM_CONDITIONAL(HAVE_AVX, test "x${have_avx}" = "xyes")
AM_CONDITIONAL(HAVE_SSE, test "x${have_sse}" = "xyes")
AM_CONDITIONAL(HAVE_PS2PDF, test "x${have_ps2pdf}" = "xyes")
AM_CONDITIONAL(FAST_BLO, test "x${use_fast_blo}" = "xyes")
AM_PROG_CC_C_O

# check for PLL
AC_CHECK_LIB([pll], [pll_partition_create], [], [AC_MSG_ERROR([PLL not found])])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/tree/Makefile
                 src/msa/Makefile
                 src/optimize/Makefile
                 src/binary/Makefile
                 src/util/Makefile
                 src/algorithm/Makefile])
AC_OUTPUT

AC_MSG_RESULT([
  $PACKAGE $VERSION

    Target:    $host_os $host_cpu
    Compiler:  ${CC}
    CFLAGS:    ${CFLAGS} ${CPPFLAGS}
    LIBS:      ${LIBS} ${LDFLAGS}

  Continue with 'make' command
])
