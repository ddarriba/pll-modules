lib_LTLIBRARIES=libpll_binary.la

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

AM_CFLAGS=-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3

LIBPLLHEADERS=\
pll.h

libpll_binary_la_SOURCES=\
     pll_binary.c \
     binary_io_operations.c \
		 ../pllmod_common.c

libpll_binary_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_binary_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_binary_la_CPPFLAGS = $(PLL_CFLAGS)
else
  libpll_binary_la_CPPFLAGS = -I$(includedir)/libpll
endif


pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll_binary.h
EXTRA_DIST = ../pllmod_common.h binary_io_operations.h
