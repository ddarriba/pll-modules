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

for mod in util msa binary tree optimize algorithm;
do

  cd $mod
  make
  cd ..

done

cd ..

cp -v --preserve libs/libpll/src/.libs/*.so* $INSTALL_DIR

cp -v --preserve src/*/*.so $INSTALL_DIR

cp -v --preserve libs/libpll/src/pll.h $INSTALL_DIR

cp -v --preserve src/*/pll_*.h src/*/pllmod_*.h $INSTALL_DIR
