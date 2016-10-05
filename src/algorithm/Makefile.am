lib_LTLIBRARIES=libpll_algorithm.la

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
AM_CPPFLAGS=-I../optimize -I../algorithm -I../tree -I../msa -I../util

LIBPLLHEADERS=\
pll.h

libpll_algorithm_la_SOURCES=\
     pllmod_algorithm.c \
     algo_callback.c \
     algo_search.c \
		 ../pllmod_common.c

libpll_algorithm_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)

# This is the version information for libtool and not the libpll version. For
#  more information see:
# https://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
libpll_algorithm_la_LDFLAGS = -version-info 0:0:0

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pllmod_algorithm.h 
EXTRA_DIST = ../pllmod_common.h algo_callback.h
