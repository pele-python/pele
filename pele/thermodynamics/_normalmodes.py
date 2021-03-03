from __future__ import print_function
import numpy as np

from pele.utils.hessian import sort_eigs

__all__ = ["normalmode_frequencies", "normalmodes", "logproduct_freq2"]


class NormalModeError(Exception):
    """raised if the normal mode results are not as expected

    this typically means that a minimum is actually a transition state,
    or a transition state is actually a higher order saddle
    """


def normalmode_frequencies(hessian, metric=None, eps=1e-4):
    """calculate (squared) normal mode frequencies

    Parameters
    ----------
    hessian: 2d array
        hessian matrix
    metric: 2d array
        mass weighted metric tensor

    Returns
    -------
    sorted array of normal mode frequencies

    """
    A = hessian
    if metric is not None:
        A = np.dot(np.linalg.pinv(metric), hessian)

    frq = np.linalg.eigvals(A)

    if np.max(np.abs(np.imag(frq))) > eps:
        print(frq)
        raise ValueError("imaginary eigenvalue in frequency calculation"
                         ", check hessian + metric tensor\n"
                         "the largest imaginary part is %g" %
                         np.max(np.abs(np.imag(frq))))

    return np.sort(np.real(frq))


def normalmodes(hessian, metric=None, eps=1e-4, symmetric=False):
    """calculate (squared) normal mode frequencies and normal mode vectors

    Parameters
    ----------
    hessian: array
        hessian marix
    metric: array
        mass weighted metric tensor
    symmetric: bool
        If true, the Hessian times the metric tensor is assumed to be symmetric.  This is
        not usually the case, even if the metric tensor is symmetric.  It is
        true if the metric tensor is the identity.

    Returns
    -------
    freq, evecs tuple array of squared frequencies and normal modes

    """
    if metric is None:
        A = hessian
        symmetric = True
    else:
        A = np.dot(np.linalg.pinv(metric), hessian)

    if symmetric:
        freq, evecs = np.linalg.eigh(A)
    else:
        freq, evecs = np.linalg.eig(A)

    if np.max(np.abs(np.imag(freq))) > eps:
        print(freq)
        raise ValueError("imaginary eigenvalue in frequency calculation"
                         ", check hessian + metric tensor\nthe largest imaginary part is %g" % np.max(
            np.abs(np.imag(freq))))

    freq = np.real(freq)
    freq, evecs = sort_eigs(freq, evecs)
    return freq, evecs


def logproduct_freq2(freqs, nzero, nnegative=0, eps=1e-4):
    """ calculate the log product of positive (squared) frequencies

    calculates
    log(product_i f_i^2)

    Parameters
    ----------
    freqs:
        array of (squared) normalmode frequencies
    nzero: int
        expected number of zero eigenvalues
    nnegative: int, optional
        expected number of negative frequencies, 0 for minimum, 1 for transition states
    eps: float, optional
        cutoff to determine if eigenvalue is no zero

    Returns
    -------
    tuple of number of considered frequencies and log product of frequencies
    """
    negative_eigs = []
    zero_eigs = []
    lnf = 0.
    n = 0
    for f in freqs:
        if np.abs(f) < eps:
            zero_eigs.append(f)
            continue
        if f < 0.:
            negative_eigs.append(f)
            continue
        lnf += np.log(f)
        n += 1

    izero = len(zero_eigs)
    inegative = len(negative_eigs)

    if nzero != izero:
        raise NormalModeError("the number of zero eigenvalues (%d) differs from the expected value (%d)" % (izero, nzero))

    if nnegative != inegative:
        if 0 < nnegative < inegative:
            raise NormalModeError("the number of negative eigenvalues (%d) is greater than expected "
                                  "(%d).  Is this a higher order saddle point?" % (inegative, nnegative))
        else:
            raise NormalModeError("the number of negative eigenvalues (%d) differs from the expected "
                                  "number (%d) (not a minimum / transition state?)" % (inegative, nnegative))

    return n, lnf

