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

OPTFLAGS=
if FAST_BLO
  OPTFLAGS+=-DNOCHECK_PERBRANCH_IMPR
endif

AM_CFLAGS=-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3

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

libpll_optimize_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS) $(OPTFLAGS)

# This is the version information for libtool and not the libpll version. For
#  more information see:
# https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
libpll_optimize_la_LDFLAGS = -version-info 0:0:0

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll_optimize.h 
EXTRA_DIST = ../pllmod_common.h lbfgsb/lbfgsb.h
