#!/bin/bash
cd source
if [ -d "build" ];then
    rm -r build
fi
cmake -B build
cd build
make
mv mymd ../..
make clean
cd ..
rm -r build
cd ..