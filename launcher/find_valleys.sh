#!/bin/bash

cd ../..

g++ ./CRT-source/find_valleys.C `root-config --cflags --glibs` -lSpectrum -o ./x/find_valleys.x

./x/find_valleys.x ./data/step2/$1_s2.root ./data/step3/$1_s3_valleys.root $1 $2
