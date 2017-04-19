#!/bin/sh

if [ "$#" -eq 0 ]; then
  PREFIX=
  PREFIX_ARG=
else
  mkdir -p $1
  PREFIX="$(cd $1 && pwd -P)" 
  PREFIX_ARG=--prefix=$PREFIX 
  PLLMOD_PREFIX_ARG="$PREFIX_ARG CPPFLAGS="-I$PREFIX/include/libpll" LDFLAGS="-L$PREFIX/lib""
fi

if [ "$#" -ge 2 ]; then
  PLL_ARGS=$2
else
  PLL_ARGS=
fi

if [ "$#" -ge 3 ]; then
  PLLMOD_ARGS=$3
else
  PLLMOD_ARGS=
fi

# configure & install libpll
cd libs/libpll
./autogen.sh
./configure $PREFIX_ARG $PLL_ARGS
make install

# configure & install modules
cd ../..
autoreconf --force --install
./configure $PLLMOD_PREFIX_ARG $PLLMOD_ARGS
make install

