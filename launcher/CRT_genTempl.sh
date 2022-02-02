#!/bin/bash

cd ../..

g++ ./CRT-source/CRT_genTemplate.C `root-config --cflags --glibs`  -lSpectrum -o ./x/CRT_genTemplate.x

./x/CRT_genTemplate.x 0 0 0 0

