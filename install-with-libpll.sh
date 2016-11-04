#!/bin/sh

if [ "$#" -eq 0 ]; then
  PREFIX=
  PREFIX_ARG=
else
  PREFIX=$(readlink -f $1)
  PREFIX_ARG=--prefix=$PREFIX 
  PLLMOD_PREFIX_ARG="$PREFIX_ARG CPPFLAGS=\"-I$PREFIX/include/libpll\" LDFLAGS=\"-L$PREFIX/lib\""
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
./configure $PLLMOD_PREFIX_ARG
make install

