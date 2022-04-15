#!/bin/bash

cd ../..

g++ ./CRT-source/optim.C `root-config --cflags --glibs`  -lSpectrum -o ./x/optim.x

for i in {0..40}
do
  export cf=$(awk "BEGIN {print 0.13 + 0.001*$i}")
  echo Constant Fraction: $cf
  ./x/optim.x 0 0 0 0 $cf
done

