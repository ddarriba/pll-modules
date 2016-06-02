#!/bin/bash

if [ "$#" -eq 0 ]; then
  INSTALL_DIR=./install
else
  INSTALL_DIR=$1
fi

if [ ! -d $INSTALL_DIR ]; then
  mkdir $INSTALL_DIR
fi


cd libs/libpll/src
make

cd ../../../src 

for mod in `ls .`;
do

  cd $mod
  make
  cd ..

done

cd ..

cp -v libs/libpll/src/*.so src/*/*.so $INSTALL_DIR

cp -v libs/libpll/src/pll.h $INSTALL_DIR

cp -v src/*/pll_*.h $INSTALL_DIR