lib_LTLIBRARIES=libpll_tree.la

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

AM_YFLAGS = -d -p `${SED} -n 's/.*_\(.*\)/pllmod_\1_/p' <<<"$*"`
AM_LFLAGS = -P `${SED} -n 's/.*_\(.*\)/pllmod_\1_/p' <<<"$*"` -o lex.yy.c

LIBPLLHEADERS=\
pll.h

libpll_tree_la_SOURCES=\
     pll_tree.c \
		 rtree_operations.c \
		 utree_operations.c \
		 utree_distances.c \
		 treeinfo.c \
		 consensus.c \
		 tree_hashtable.c \
		 split_utree.y \
		 lex_split.l \
		 ../pllmod_common.c

libpll_tree_la_CFLAGS = $(AM_CFLAGS) $(AVXFLAGS) $(SSEFLAGS)
libpll_tree_la_LDFLAGS = -version-info 0:0:0
if HAVE_PLL_DPKG
  libpll_tree_la_CPPFLAGS = $(PLL_CFLAGS)
else
  libpll_tree_la_CPPFLAGS = -I$(includedir)/libpll
endif

pkgincludedir=$(includedir)/libpll
pkginclude_HEADERS = pll_tree.h ../pllmod_common.h
EXTRA_DIST = ../pllmod_common.h tree_hashtable.h
