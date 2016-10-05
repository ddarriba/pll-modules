#!/bin/sh

PLLMOD_ROOT=`dirname $(readlink -f $0)`
PLL_ROOT=$PLLMOD_ROOT/libs/libpll

if [ "$#" -eq 0 ]; then
  PREFIX=
  PREFIX_ARG=
else
  PREFIX=$(readlink -f $1)
  PREFIX_ARG=--prefix=$PREFIX
  mkdir -p $PREFIX
fi

# configure & install libpll
cd libs/libpll
./autogen.sh
./configure $PREFIX_ARG
make install

# configure & install modules
cd ../..
autoreconf --force --install
./configure $PREFIX_ARG CPPFLAGS="-I$PLL_ROOT/src -L$PLL_ROOT/src/.libs" LDFLAGS="-L$PLL_ROOT/src/.libs"
make install

