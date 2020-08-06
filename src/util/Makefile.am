lib_LTLIBRARIES=libpll_util.la

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

libpll_util_la_SOURCES=\
     models.c \
     models_dna.c \
     models_aa.c \
     models_gt.c \
     models_mult.c \
		 ../pllmod_common.c

libpll_util_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_util_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_util_la_CPPFLAGS = $(PLL_CFLAGS)
else
  libpll_util_la_CPPFLAGS = -I$(includedir)/libpll
endif

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pllmod_util.h ../pllmod_common.h
EXTRA_DIST = ../pllmod_common.h
