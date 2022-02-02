#!/bin/bash

cd ../..

g++ ./CRT-source/CRT_genCutsG.C `root-config --cflags --glibs` -lSpectrum -o ./x/CRT_genCutsG.x

./x/CRT_genCutsG.x ./data/step2/$1_s2.root ./data/step3/$1_cutg.root $1 $2
