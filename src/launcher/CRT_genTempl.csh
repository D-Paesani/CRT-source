#!/bin/tcsh

cd ../..

g++ ./src/CRT_genTemplate.C `root-config --cflags --glibs` -o ./src/x/CRT_genTemplate.x

./src/x/CRT_genTemplate.x 0 0 0 0

