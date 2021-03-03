from __future__ import print_function
import unittest
import numpy as np

from pele.systems import LJCluster
from pele.optimize import MYLBFGS, LBFGS, lbfgs_py

class TestMYLBFGS_LBFGS(unittest.TestCase):
    def setUp(self):
        self.setUp1(verbose=False)

    def setUp1(self, verbose=False, **kwargs):
        np.random.seed(0)
        natoms = 18
        self.system = LJCluster(natoms)
        self.pot = self.system.get_potential()
        x = self.system.get_random_configuration()
        ret = lbfgs_py(x, self.pot, tol=10)
        self.x = ret.coords
        
        self.kwargs = kwargs
        self.verbose = verbose
    
        self.M = 4
        if self.verbose: iprint=1
        else: iprint = -1
        self.myo = MYLBFGS(self.x, self.pot, iprint=iprint, debug=True, M=self.M)
        self.o = LBFGS(self.x, self.pot, iprint=iprint, debug=True, M=self.M, **self.kwargs)

    def test(self):
        N = self.x.size
        M = self.M
        myo = self.myo
        o = self.o

        # do one iteration
        for i in range(3 * self.M):   
            myo.one_iteration()
            o.one_iteration()
            if self.verbose:
                print("")
                print("H0", myo.H0, o.H0)
                print("rho  ", o.rho[:])
                print("myrho", myo.W[N:N+M])

        myret = myo.get_result()
        ret = o.get_result()
        
        self.assertAlmostEqual(ret.energy, myret.energy, 4)
        self.assertLess(np.max(np.abs(myret.coords - ret.coords)), 1e-6)
    
        # do a second iteration
        for i in range(1):
            myo.one_iteration()
            o.one_iteration()
        myret = myo.get_result()
        ret = o.get_result()
        
        if self.verbose:
            print("H0", myret.H0, ret.H0)
            print("rho  ", o.rho[:])
            print("myrho", myo.W[N:N+M])

        self.assertAlmostEqual(ret.energy, myret.energy, 4)
        self.assertLess(np.max(np.abs(myret.coords - ret.coords)), 1e-6)
    
    def test_complete(self):
        myret = self.myo.run()
        ret = self.o.run()
        self.assertEqual(ret.nfev, myret.nfev)
        self.assertEqual(ret.nsteps, myret.nsteps)
        self.assertAlmostEqual(ret.energy, myret.energy, 4)
        self.assertLess(np.max(np.abs(myret.coords - ret.coords)), 1e-6)
#        print "ediff", ret.energy - myret.energy
        

class TestMYLBFGS_LBFGS_Cython(TestMYLBFGS_LBFGS):
    def setUp1(self, verbose=False, **kwargs):
        TestMYLBFGS_LBFGS.setUp1(self, verbose=verbose, **kwargs)
        self.o._cython=True

class TestMYLBFGS_LBFGS_fortran(TestMYLBFGS_LBFGS):
    def setUp(self):
        self.setUp1(fortran=True)

if __name__ == "__main__":
    unittest.main()

        
