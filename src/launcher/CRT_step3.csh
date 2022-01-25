#!/bin/tcsh

cd ../..

g++ ./src/CRT_step3.C `root-config --cflags --glibs` -o ./src/x/CRT_step3.x

./src/x/./CRT_step3.x ./data/step2/$1_s2.root ./data/step3/$1_s3.root $1 $2

