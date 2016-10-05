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
		 ../pllmod_common.c

libpll_util_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)

# This is the version information for libtool and not the libpll version. For
#  more information see:
# https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
libpll_util_la_LDFLAGS = -version-info 0:0:0

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pllmod_util.h
EXTRA_DIST = ../pllmod_common.h
