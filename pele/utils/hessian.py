"""
Tools for manipulating the Hessian.  In particular, for finding eigenvalues and eigenvectors

.. currentmodule:: pele.utils.hessian

.. autosummary::
    :toctree: generated
    
    get_eig
    get_eigvals
    get_sorted_eig
    get_smallest_eig
    make_sparse

"""
import numpy as np 
from pele.potentials.lj import LJ

__all__ = ["get_eig", "get_eigvals", "get_sorted_eig", "get_smallest_eig", "make_sparse"]

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

def sort_eigs(evals, evecs):
    """return the sorted eigenvalues and eigenvectors"""
    mylist = [(evals[i], evecs[:,i]) for i in range(len(evals))]
    sortlist = sorted(mylist, key=lambda x:x[0])
    evals = np.array([wv[0] for wv in sortlist])
    evecs = evecs.copy()
    for i in range(len(evals)):
        evecs[:,i] = sortlist[i][1]
    return evals, evecs
    

def get_sorted_eig(hess, **kwargs):
    """return the sorted eigenvalues and eigenvectors of a Hessian sorted"""
    evals, evecs = get_eig(hess, **kwargs)
    # now sort them
    try:
        mylist = [(evals[i], evecs[:,i]) for i in range(len(evals))]
        sortlist = sorted(mylist, key=lambda x:x[0])
    except ValueError:
        import sys
        print >> sys.stderr, "evals, evecs", evals.shape, evecs.shape
        print >> sys.stderr, "evals", evals
        print >> sys.stderr, "evecs", evecs
        print >> sys.stderr, evals[0], evecs[:,1]
        print >> sys.stderr, mylist
        raise
    evals = np.array([wv[0] for wv in sortlist])
    for i in range(len(evals)):
        evecs[:,i] = sortlist[i][1]
    return evals, evecs

def get_smallest_eig(hess, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian"""
    evals, evecs = get_sorted_eig(hess, **kwargs)
    return evals[0], evecs[:,0].flatten()

def get_smallest_eig_arpack(hess, tol=1e-3, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian
    
    use arpack 
    """
    import scipy.sparse
    from scipy.sparse.linalg import eigsh, eigs
    from scipy.sparse.linalg.eigen.arpack.arpack import ArpackNoConvergence
    import sys
    try:
        e, v = eigsh(hess, which="SA", k=1, maxiter=1000, tol=tol)
    except ArpackNoConvergence as err:
        sys.stderr.write("ArpackNoConvergence raised\n")
        if scipy.sparse.issparse(hess):
            hess = hess.todense()
        return get_smallest_eig(hess, **kwargs)
    return e[0], v[:,0].flatten()

def get_smallest_eig_sparse(hess, cutoff=1e-1, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian
    
    use arpack, and set all hessian values less than cutoff to zero
    """
    import scipy.sparse.linalg
    newhess = np.where(np.abs(hess) < cutoff, 0., hess)
    # i can't get it to work taking only the upper or lower triangular matrices
#    sparsehess = scipy.sparse.tril(newhess, format="csr")
    sparsehess = scipy.sparse.csr_matrix(newhess)
#    print "dense len", len(hess.reshape(-1)), "sparse len", len(sparsehess.nonzero()[0])
    return get_smallest_eig_arpack(sparsehess, **kwargs)

def get_smallest_eig_nohess(coords, system, **kwargs):
    """find the smallest eigenvalue and eigenvector without a hessian
    
    this is just a wrapper for findLowestEigenVector
    
    See Also
    --------
    pele.transition_states.findLowestEigenVector
    """
    from pele.transition_states import findLowestEigenVector
    ret = findLowestEigenVector(coords, system.get_potential(), orthogZeroEigs=system.get_orthogonalize_to_zero_eigenvectors(), **kwargs)
    return ret.eigenval, ret.eigenvec

def make_sparse(hess, **kwargs):
    """return a sparse form of the hessian using scipy.sparse
    
    this function returns the hessian in compressed sparse column (CSC) form
    
    Advantages of the CSC format:
        - efficient arithmetic operations CSC + CSC, CSC * CSC, etc.
        - efficient column slicing
        - fast matrix vector products (CSR, BSR may be faster)
    Disadvantages of the CSC format:
        - slow row slicing operations (consider CSR)
        - changes to the sparsity structure are expensive (consider LIL or DOK)
    """
    import scipy.sparse as sparse
    return sparse.csc_matrix(hess)

#
# only testing stuff below here
#    
    
    
import unittest
class TestEig(unittest.TestCase):
    def setUp(self):
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
        w, v = get_smallest_eig_nohess(self.x, self.system, tol=1e-9, dx=1e-3)
#        print vs.shape, v.shape
        self.assertAlmostEqual(ws, w, 1)
        dot = np.dot(v, vs) / (np.linalg.norm(v) * np.linalg.norm(vs))
        dot = np.abs(dot)
        self.assertAlmostEqual(dot, 1., 1)

        

def size_scaling_smallest_eig(natoms):
    from pele.systems import LJCluster
    import time, sys
    system = LJCluster(natoms)
    pot = system.get_potential()
    quencher = system.get_minimizer(tol=10.)
    
    time1 = 0.
    time2 = 0.
    time3 = 0.
    time4 = 0.
    for i in range(100):
        coords = system.get_random_configuration()
#        print "len(coords)", len(coords)
        coords = quencher(coords)[0]
        e, g, h = pot.getEnergyGradientHessian(coords)
        
        t0 = time.time()
        w1, v1 = get_smallest_eig(h)
        t1 = time.time()
        w, v = get_smallest_eig_arpack(h)
        t2 = time.time()
        w2, v2 = get_smallest_eig_sparse(h)
        t3 = time.time()
        w3, v3 = get_smallest_eig_nohess(coords, system, tol=1e-3)
        t4 = time.time()
        
        time1 += t1-t0
        time2 += t2-t1
        time3 += t3-t2
        time4 += t4-t3
        
        wdiff = np.abs(w-w1) / np.max(np.abs([w,w1]))
        if wdiff > 5e-3:
            sys.stderr.write("eigenvalues for dense  are different %g %g normalized diff %g\n" % (w1, w, wdiff))
        wdiff = np.abs(w-w2) / np.max(np.abs([w,w2]))
        if wdiff > 5e-2:
            sys.stderr.write("eigenvalues for sparse are different %g %g normalized diff %g\n" % (w2, w, wdiff))
        wdiff = np.abs(w-w3) / np.max(np.abs([w,w3]))
        if wdiff > 5e-2:
            sys.stderr.write("eigenvalues for nohess are different %g %g normalized diff %g\n" % (w3, w, wdiff))
#    print "times", n, t1-t0, t2-t1, w1, w
    print "times", n, time1, time2, time3, time4
    sys.stdout.flush()

def plot_hist(hess):
    import pylab as pl
    pl.hist(np.log10(np.abs(hess.reshape(-1))))
    pl.show()

if __name__ == "__main__":
    from pele.systems import LJCluster
    natoms = 30
    system = LJCluster(natoms)
    pot = system.get_potential()
    coords = system.get_random_configuration()
    
    xmin = system.get_random_minimized_configuration()[0]
    e, g, h = pot.getEnergyGradientHessian(xmin)
    evals = get_eigvals(h)
    print evals

    quencher = system.get_minimizer(tol=10.)
    coords = quencher(coords)[0]
    e, g, h = pot.getEnergyGradientHessian(coords)
    w1, v1 = get_smallest_eig(h)
    print w1
    w, v = get_smallest_eig_arpack(h)
    print w
    w2, v2 = get_smallest_eig_sparse(h)
    print w2, w2/w1
    w3, v3 = get_smallest_eig_nohess(coords, system)
    print w3, w3/w1
#    plot_hist(h)
#    exit()
    
    if False:
        n = 10
        while n < 500:
            size_scaling_smallest_eig(int(n))
            n *= 1.2
    
    unittest.main()
    
