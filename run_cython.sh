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
echo "Compiling _harmonic_cpp.pyx"
cython --cplus pele/potentials/_harmonic_cpp.pyx 
echo "Compiling _modified_fire_cpp.pyx"
cython --cplus -a pele/optimize/_modified_fire_cpp.pyx
echo "Compiling _pele_mc.pyx"
cython --cplus -a playground/monte_carlo/_pele_mc.pyx
echo "Compiling _monte_carlo_cpp.pyx"
cython --cplus -a playground/monte_carlo/_monte_carlo_cpp.pyx
echo "Compiling _takestep_cpp.pyx"
cython --cplus -a playground/monte_carlo/_takestep_cpp.pyx
echo "Compiling _accept_test_cpp.pyx"
cython --cplus -a playground/monte_carlo/_accept_test_cpp.pyx
echo "Compiling _conf_test_cpp.pyx"
cython --cplus -a playground/monte_carlo/_conf_test_cpp.pyx
echo "Compiling _action_cpp.pyx"
cython --cplus -a playground/monte_carlo/_action_cpp.pyx