#!/bin/bash -e

echo "Compiling _pele.pyx"
cython --cplus pele/potentials/_pele.pyx 
echo "Compiling _lj.pyx"
cython --cplus pele/potentials/_lj_cpp.pyx 
#echo "Compiling _lj_cython.pyx"
#cython --cplus _lj_cython.pyx 
echo "Compiling _lbfgs_cpp.pyx"
cython --cplus -a pele/optimize/_lbfgs_cpp.pyx 
echo "Compiling _pythonpotential.pyx"
cython --cplus -a pele/potentials/_pythonpotential.pyx 
echo "Compiling _morse_cpp.pyx"
cython --cplus pele/potentials/_morse_cpp.pyx 
echo "Compiling _hs_wca_cpp.pyx"
cython --cplus pele/potentials/_hs_wca_cpp.pyx 
echo "Compiling _wca_cpp.pyx"
cython --cplus pele/potentials/_wca_cpp.pyx 
echo "Compiling _modified_fire_cpp.pyx"
cython --cplus -a pele/optimize/_modified_fire_cpp.pyx 