#!/bin/bash

rm benchmark.dat

for natoms in 2 5 10 13 20 30 38 40 72 100; do
  echo $natoms
  python test.py $natoms 100000 | tail -1 >> benchmark.dat
done
