import unittest
import numpy as np

from pele.utils.hessian import *
from pele.utils.hessian import get_smallest_eig_nohess, get_smallest_eig_sparse, get_smallest_eig_arpack

class TestEig(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        from pele.systems import LJCluster
        natoms = 10
        self.system = LJCluster(natoms)
        system = self.system
        self.pot = system.get_potential()
        quencher = system.get_minimizer(tol=2.)
        x = system.get_random_configuration()
        ret = quencher(x)
        self.x = ret[0]

        self.xmin = system.get_random_minimized_configuration()[0]
        
        e, g, self.h = self.pot.getEnergyGradientHessian(self.x)
        e, g, self.hmin = self.pot.getEnergyGradientHessian(self.xmin)
        
    
    def numerical_eig_from_vec(self, x, vec, eps=1e-6):
        x = x.copy()
        x += vec * eps
        eplus, gplus = self.pot.getEnergyGradient(x)
        x -= 2. * vec * eps
        eminus, gminus = self.pot.getEnergyGradient(x)
        eval = np.dot((gplus - gminus), vec) / (2. * eps)
        return eval
    
    def test_minimum(self):
        w = get_eigvals(self.hmin)
        wmin = np.min(w)
        self.assertGreater(wmin, -1e-5)

    def test_eig_eigval(self):
        w0 = get_eigvals(self.h)
        w, v = get_eig(self.h)
        diff = np.max(np.abs(w-w0))
        self.assertLess(diff, 1e-5)

    def test_numeric(self):
        wlist, vlist = get_eig(self.h)
        eps = 1e-6
        for i in range(len(wlist)):
            w = wlist[i]
            v = vlist[:,i]
            eval = self.numerical_eig_from_vec(self.x, v)
            self.assertAlmostEqual(w, eval, 5)

    def test_numeric_sorted(self):
        wlist, vlist = get_sorted_eig(self.h)
        eps = 1e-6
        for i in range(len(wlist)):
            w = wlist[i]
            v = vlist[:,i]
            eval = self.numerical_eig_from_vec(self.x, v)
            self.assertAlmostEqual(w, eval, 5)
            
    def test_sorting(self):
        w, v = get_eig(self.h)
        ws, vs = get_sorted_eig(self.h)
        wsort = np.array(sorted(w))
        diff = np.max(np.abs(ws - wsort))
        self.assertLess(diff, 1e-5)
        
#        print "unsorted", v
#        print "sorted", vs
        isort = sorted([(w[i], i) for i in range(len(w))])
        indices = [i for wval, i in isort]
        for i, j in enumerate(indices):
            self.assertAlmostEqual(ws[i], w[j], 5)
            if np.abs(w[j]) > .01:
#                print w[j]
                v1 = vs[:,i]
                v2 = v[:,j]
                dot = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                self.assertAlmostEqual(dot, 1., 5)
                diff = np.max(np.abs(vs[:,i] - v[:,j]))
                self.assertLess(diff, 1e-5)

    def test_smallest_eig(self):
        ws, vs = get_sorted_eig(self.h)
        ws = ws[0]
        vs = vs[:,0]
        w, v = get_smallest_eig(self.h)
        self.assertAlmostEqual(ws, w, 6)
        dot = np.dot(v, vs) / (np.linalg.norm(v) * np.linalg.norm(vs))
        self.assertAlmostEqual(dot, 1., 5)

    def test_smallest_eig1(self):
        ws, vs = get_smallest_eig(self.h)
        w, v = get_smallest_eig_arpack(self.h, tol=1e-9)
        self.assertAlmostEqual(ws, w, 3)
        dot = np.dot(v, vs) / (np.linalg.norm(v) * np.linalg.norm(vs))
        dot = np.abs(dot)
        self.assertAlmostEqual(dot, 1., 3)

    def test_smallest_eig2(self):
        ws, vs = get_smallest_eig(self.h)
        w, v = get_smallest_eig_sparse(self.h, cutoff=1e-2, tol=1e-9)
#        print vs.shape, v.shape
        self.assertAlmostEqual(ws, w, 2)
        dot = np.dot(v, vs) / (np.linalg.norm(v) * np.linalg.norm(vs))
        dot = np.abs(dot)
        self.assertAlmostEqual(dot, 1., 2)

    def test_smallest_eig_nohess(self):
        ws, vs = get_smallest_eig(self.h)
        w, v = get_smallest_eig_nohess(self.x, self.system, tol=1e-9, dx=1e-6)
#        print vs.shape, v.shape
        self.assertAlmostEqual(ws, w, 1)
        dot = np.dot(v, vs) / (np.linalg.norm(v) * np.linalg.norm(vs))
        dot = np.abs(dot)
        self.assertAlmostEqual(dot, 1., 1)
        
if __name__ == "__main__":
    unittest.main()

