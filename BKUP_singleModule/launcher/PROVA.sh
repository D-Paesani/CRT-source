#!/bin/bash

cd ../..

g++ ./CRT-source/PROVA.C `root-config --cflags --glibs`  -lSpectrum -o ./x/PROVA.x

./x/PROVA.x 0 0 0 0

