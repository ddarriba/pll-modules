lib_LTLIBRARIES=libpll_msa.la

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

libpll_msa_la_SOURCES=\
     pll_msa.c \
		 ../pllmod_common.c

libpll_msa_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_msa_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_msa_la_CPPFLAGS = $(PLL_CFLAGS)
else
  libpll_msa_la_CPPFLAGS = -I$(includedir)/libpll
endif


pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll_msa.h 
EXTRA_DIST = ../pllmod_common.h
