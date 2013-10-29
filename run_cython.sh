#!/bin/bash -e

echo "Compiling _pele.pyx"
cython --cplus pele/potentials/_pele.pyx 
echo "Compiling _lj.pyx"
cython --cplus pele/potentials/_lj_cpp.pyx 
#echo "Compiling _lj_cython.pyx"
#cython --cplus _lj_cython.pyx 
echo "Compiling _lbfgs.pyx"
cython --cplus -a pele/optimize/_lbfgs_cpp.pyx 
echo "Compiling _pythonpotential.pyx"
cython --cplus -a pele/potentials/_pythonpotential.pyx 

