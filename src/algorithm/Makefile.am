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

AM_CFLAGS=-Wall -Wsign-compare -D_GNU_SOURCE -std=c99 -O3 #-DDEBUG
AM_CPPFLAGS=-I../optimize -I../algorithm -I../tree -I../msa -I../util

LIBPLLHEADERS=\
pll.h

libpll_algorithm_la_SOURCES=\
     pllmod_algorithm.c \
     algo_callback.c \
     algo_search.c \
		 ../pllmod_common.c

libpll_algorithm_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_algorithm_la_CPPFLAGS = $(AM_CPPFLAGS) 
libpll_algorithm_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_algorithm_la_CPPFLAGS += $(PLL_CFLAGS)
else
  libpll_algorithm_la_CPPFLAGS += -I$(includedir)/libpll
endif


pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pllmod_algorithm.h 
EXTRA_DIST = ../pllmod_common.h algo_callback.h
