"""
Tools for manipulating the Hessian.  In particular, for finding eigenvalues and eigenvectors

.. currentmodule:: pygmin.utils.hessian

.. autosummary::
    :toctree: generated
    
    get_eig
    get_eigvals
    get_sorted_eig
    get_smallest_eig

"""
import numpy as np 
from pygmin.potentials.lj import LJ

__all__ = ["get_eig", "get_eigvals", "get_sorted_eig", "get_smallest_eig"]

def get_eigvals(hess, **kwargs):
    """return the eigenvalues of a Hessian (symmetric)
    
    The following docs are from numpy.linalg.eigvalsh
    
    Compute the eigenvalues of a Hermitian or real symmetric matrix.

    Main difference from eigh: the eigenvectors are not computed.

    Parameters
    ----------
    a : array_like, shape (M, M)
        A complex- or real-valued matrix whose eigenvalues are to be
        computed.
    UPLO : {'L', 'U'}, optional
        Specifies whether the calculation is done with the lower triangular
        part of `a` ('L', default) or the upper triangular part ('U').

    Returns
    -------
    w : ndarray, shape (M,)
        The eigenvalues, not necessarily ordered, each repeated according to
        its multiplicity.

    Raises
    ------
    LinAlgError
        If the eigenvalue computation does not converge.

    """ 
    return np.linalg.eigvalsh(hess, **kwargs)


def get_eig(hess, **kwargs):
    """return the eigenvalue and eigenvectors of a Hessian (symmetric)
    
    The following is from numpy.linalg.eigh
    
    Return the eigenvalues and eigenvectors of a Hermitian or symmetric matrix.

    Returns two objects, a 1-D array containing the eigenvalues of `a`, and
    a 2-D square array or matrix (depending on the input type) of the
    corresponding eigenvectors (in columns).

    Parameters
    ----------
    a : array_like, shape (M, M)
        A complex Hermitian or real symmetric matrix.
    UPLO : {'L', 'U'}, optional
        Specifies whether the calculation is done with the lower triangular
        part of `a` ('L', default) or the upper triangular part ('U').

    Returns
    -------
    w : ndarray, shape (M,)
        The eigenvalues, not necessarily ordered.
    v : ndarray, or matrix object if `a` is, shape (M, M)
        The column ``v[:, i]`` is the normalized eigenvector corresponding
        to the eigenvalue ``w[i]``.

    Raises
    ------
    LinAlgError
        If the eigenvalue computation does not converge.

    """ 
    return np.linalg.eigh(hess, **kwargs)

def get_sorted_eig(hess, **kwargs):
    """return the sorted eigenvalues and eigenvectors of a Hessian"""
    evals, evecs = get_eig(hess, **kwargs)
    # now sort them
    sortlist = sorted([(evals[i], evecs[:,i]) for i in range(len(evals))])
    
    evals = np.array([wv[0] for wv in sortlist])
    for i in range(len(evals)):
        evecs[:,i] = sortlist[i][1]
    return evals, evecs

def get_smallest_eig(hess, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian"""
    evals, evecs = get_smallest_eig(hess, **kwargs)
    return eval[0], evecs[:,0]
    
    
import unittest
class TestEig(unittest.TestCase):
    def setUp(self):
        from pygmin.systems import LJCluster
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


if __name__ == "__main__":
#    from pygmin.systems import LJCluster
#    natoms = 10
#    system = LJCluster(10)
#    pot = system.get_potential()
#    coords = system.get_random_configuration()
#    
#    xmin = system.get_random_minimized_configuration()[0]
#    e, g, h = pot.getEnergyGradientHessian(xmin)
#    evals = get_eigvals(h)
#    print evals
    
    unittest.main()
    