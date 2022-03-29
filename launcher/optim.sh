#!/bin/bash

cd ../..

g++ ./CRT-source/optim.C `root-config --cflags --glibs`  -lSpectrum -o ./x/optim.x

./x/optim.x 0 0 0 0

