#!/bin/bash

echo "Compiling _pygmin.pyx"
cython --cplus _pygmin.pyx 
echo "Compiling _lj.pyx"
cython --cplus _lj.pyx 
echo "Compiling _lj_cython.pyx"
cython --cplus _lj_cython.pyx 

