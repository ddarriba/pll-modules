lib_LTLIBRARIES=libpll_optimize.la

if HAVE_AVX
AVXFLAGS=-mavx
else
AVXFLAGS=
endif

if HAVE_SSE
SSEFLAGS=-msse3
else
SSEFLAGS=
endif

AM_CFLAGS=-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3 #-DDEBUG

LIBPLLHEADERS=\
pll.h

libpll_optimize_la_SOURCES=\
     pll_optimize.c \
     opt_algorithms.c \
     lbfgsb/lbfgsb.c \
     lbfgsb/linesearch.c \
     lbfgsb/linpack.c \
     lbfgsb/miniCBLAS.c \
     lbfgsb/subalgorithms.c \
		 ../pllmod_common.c

libpll_optimize_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_optimize_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_optimize_la_CPPFLAGS = $(PLL_CFLAGS)
else
  libpll_optimize_la_CPPFLAGS = -I$(includedir)/libpll
endif


pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll_optimize.h ../pllmod_common.h
EXTRA_DIST = ../pllmod_common.h lbfgsb/lbfgsb.h
