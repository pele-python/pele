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
from __future__ import print_function
import numpy as np

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


def sort_eigs(evals, evecs, reverse=False):
    """return the sorted eigenvalues and eigenvectors"""
    mylist = [(evals[i], evecs[:, i]) for i in range(len(evals))]
    sortlist = sorted(mylist, key=lambda x: x[0], reverse=reverse)
    evals = np.array([wv[0] for wv in sortlist])
    evecs = evecs.copy()
    for i in range(len(evals)):
        evecs[:, i] = sortlist[i][1]
    return evals, evecs


def get_sorted_eig(hess, **kwargs):
    """return the sorted eigenvalues and eigenvectors of a Hessian sorted"""
    evals, evecs = get_eig(hess, **kwargs)
    # now sort them
    try:
        mylist = [(evals[i], evecs[:, i]) for i in range(len(evals))]
        sortlist = sorted(mylist, key=lambda x: x[0])
    except ValueError:
        import sys

        print("evals, evecs", evals.shape, evecs.shape, file=sys.stderr)
        print("evals", evals, file=sys.stderr)
        print("evecs", evecs, file=sys.stderr)
        print(evals[0], evecs[:, 1], file=sys.stderr)
        print(mylist, file=sys.stderr)
        raise
    evals = np.array([wv[0] for wv in sortlist])
    for i in range(len(evals)):
        evecs[:, i] = sortlist[i][1]
    return evals, evecs


def get_smallest_eig(hess, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian"""
    evals, evecs = get_sorted_eig(hess, **kwargs)
    return evals[0], evecs[:, 0].flatten()


def get_smallest_eig_arpack(hess, tol=1e-3, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian
    
    use arpack 
    """
    import scipy.sparse
    from scipy.sparse.linalg import eigsh
    from scipy.sparse.linalg.eigen.arpack.arpack import ArpackNoConvergence
    import sys

    try:
        e, v = eigsh(hess, which="SA", k=1, maxiter=1000, tol=tol)
    except ArpackNoConvergence:
        sys.stderr.write("ArpackNoConvergence raised\n")
        if scipy.sparse.issparse(hess):
            hess = hess.todense()
        return get_smallest_eig(hess, **kwargs)
    return e[0], v[:, 0].flatten()


def get_smallest_eig_sparse(hess, cutoff=1e-1, **kwargs):
    """return the smallest eigenvalue and associated eigenvector of a Hessian
    
    use arpack, and set all hessian values less than cutoff to zero
    """
    import scipy.sparse.linalg

    newhess = np.where(np.abs(hess) < cutoff, 0., hess)
    # i can't get it to work taking only the upper or lower triangular matrices
    # sparsehess = scipy.sparse.tril(newhess, format="csr")
    sparsehess = scipy.sparse.csr_matrix(newhess)
    # print "dense len", len(hess.reshape(-1)), "sparse len", len(sparsehess.nonzero()[0])
    return get_smallest_eig_arpack(sparsehess, **kwargs)


def get_smallest_eig_nohess(coords, system, **kwargs):
    """find the smallest eigenvalue and eigenvector without a hessian
    
    this is just a wrapper for findLowestEigenVector
    
    See Also
    --------
    pele.transition_states.findLowestEigenVector
    """
    from pele.transition_states import findLowestEigenVector

    ret = findLowestEigenVector(coords, system.get_potential(),
                                orthogZeroEigs=system.get_orthogonalize_to_zero_eigenvectors(), **kwargs)
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






def size_scaling_smallest_eig(natoms):  # pragma: no cover
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
        # print "len(coords)", len(coords)
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

        time1 += t1 - t0
        time2 += t2 - t1
        time3 += t3 - t2
        time4 += t4 - t3

        wdiff = np.abs(w - w1) / np.max(np.abs([w, w1]))
        if wdiff > 5e-3:
            sys.stderr.write("eigenvalues for dense  are different %g %g normalized diff %g\n" % (w1, w, wdiff))
        wdiff = np.abs(w - w2) / np.max(np.abs([w, w2]))
        if wdiff > 5e-2:
            sys.stderr.write("eigenvalues for sparse are different %g %g normalized diff %g\n" % (w2, w, wdiff))
        wdiff = np.abs(w - w3) / np.max(np.abs([w, w3]))
        if wdiff > 5e-2:
            sys.stderr.write("eigenvalues for nohess are different %g %g normalized diff %g\n" % (w3, w, wdiff))
        # print "times", n, t1-t0, t2-t1, w1, w
    print("times", n, time1, time2, time3, time4)
    sys.stdout.flush()


def plot_hist(hess):  # pragma: no cover
    import pylab as pl

    pl.hist(np.log10(np.abs(hess.reshape(-1))))
    pl.show()


def test():  # pragma: no cover
    from pele.systems import LJCluster

    natoms = 30
    system = LJCluster(natoms)
    pot = system.get_potential()
    coords = system.get_random_configuration()

    xmin = system.get_random_minimized_configuration()[0]
    e, g, h = pot.getEnergyGradientHessian(xmin)
    evals = get_eigvals(h)
    print(evals)

    quencher = system.get_minimizer(tol=10.)
    coords = quencher(coords)[0]
    e, g, h = pot.getEnergyGradientHessian(coords)
    w1, v1 = get_smallest_eig(h)
    print(w1)
    w, v = get_smallest_eig_arpack(h)
    print(w)
    w2, v2 = get_smallest_eig_sparse(h)
    print(w2, w2 / w1)
    w3, v3 = get_smallest_eig_nohess(coords, system)
    print(w3, w3 / w1)
    # plot_hist(h)
    # exit()

    if False:
        n = 10
        while n < 500:
            size_scaling_smallest_eig(int(n))
            n *= 1.2


if __name__ == "__main__":
    test()

