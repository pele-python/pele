import time
import sys
import tempfile
import subprocess
import argparse

import numpy as np

from pele.potentials._lj_cpp import LJ as LJcpp
from pele.potentials.lj import LJ as LJfortran

from pele.optimize._quench import mylbfgs, lbfgs_py, lbfgs_cpp, lbfgs_scipy

optimizers = {
              mylbfgs : "fortran lbfgs",
              lbfgs_py : "python lbfgs",
              lbfgs_cpp : "cpp lbfgs",
              lbfgs_scipy : "scipy lbfgs",
              }

potentials = {
              LJcpp : "cpp LJ",
              LJfortran : "fortran LJ"
              }

def bench_pure_cpp(prog, coords):
    print ""
    coords_file = "coords.bench"
    print "calling pure c++ executable {} using subprocess module".format(prog)
    np.savetxt(coords_file, coords.ravel())
    
    t0 = time.time()
    p = subprocess.call([prog, coords_file])
    t1 = time.time()
    
    print "pure c++ LJ minimization                : time {}".format(t1-t0)

def bench_python_cpp(coords):
    print ""
    lj = LJcpp()
    t0 = time.time()
    lbfgs_cpp(coords.copy(), lj, iprint=1000, nsteps=100000)
    t1 = time.time()
    
    print "python wrapped pure c++ LJ minimization : time {}".format(t1-t0)

    

def do_benchmarks(natoms=1000, prog="./benchmarks/bench_minimization"):
    coords = np.random.uniform(-1, 1, 3*natoms) * 1.5 * natoms**(1./3)
    lj = LJcpp()
    print "intial energy :", lj.getEnergy(coords.copy())
    
    bench_pure_cpp(prog, coords)
    bench_python_cpp(coords)
    
    opt_kwargs = dict(nsteps=100000)
    print "\n"
    for pot, pot_name in potentials.iteritems():
        for opt, opt_name in optimizers.iteritems():
            p = pot()
            t0 = time.time()
            opt(coords.copy(), p, **opt_kwargs)
            t1 = time.time()
            print "{pot_name:10}: {opt_name:13}: time {time}".format(pot_name=pot_name,
                                                               opt_name=opt_name,
                                                               time=t1-t0
                                                               )
            sys.stdout.flush()

def main():
    # guess the location of the cpp executable
    # if run from the build directory
    prog = "./benchmarks/bench_minimization"
    
    parser = argparse.ArgumentParser(description="LJ minimization benchmarks")
    parser.add_argument("--cpp_exec", help="the cpp executable", type=str, default=prog)
    parser.add_argument("--natoms", type=int, default=200)
    args = parser.parse_args()
    
    do_benchmarks(natoms=args.natoms, prog=args.cpp_exec)
    
if __name__ == "__main__":
    main() 