import numpy as np

from pygmin.utils.hessian import sort_eigs

__all__ = ["normalmode_frequencies", "normalmodes", "logproduct_freq2"]

def normalmode_frequencies(hessian, metric=None, eps=1e-4):
    '''calculate (squared) normal mode frequencies
    
    Parameters
    ----------
    hessian: 2d array
        hessian matrix
    metric: 2d array
        mass weighted metric tensor
        
    Returns
    -------
    sorted array of normal mode frequencies
    
    '''
    A = hessian
    if metric is not None:
        A = np.dot(np.linalg.pinv(metric), hessian)
   
    frq = np.linalg.eigvals(A)
    
    if(np.max(np.abs(np.imag(frq))) > eps):
        print frq
        raise ValueError("imaginary eigenvalue in frequency calculation"
                         ", check hessian + metric tensor\nthe largest imaginary part is %g"%np.max(np.abs(np.imag(frq))))
    
    return np.sort(np.real(frq))

def normalmodes(hessian, metric=None, eps=1e-4):
    '''calculate (squared) normal mode frequencies and normal mode vectors
    
    Parameters
    ----------
    hessian:
        hessian marix
    metric: 
        mass weighted metric tensor
        
    Returns
    -------
    freq, evecs tuple array of squared frequencies and normal modes
    
    '''
    A = hessian
    if metric is not None:
        A = np.dot(np.linalg.pinv(metric), hessian)
   
    freq, evecs = np.linalg.eig(A)
    
    if(np.max(np.abs(np.imag(freq))) > eps):
        print freq
        raise ValueError("imaginary eigenvalue in frequency calculation"
                         ", check hessian + metric tensor\nthe largest imaginary part is %g"%np.max(np.abs(np.imag(freq))))
    
    freq = np.real(freq)
    freq, evecs = sort_eigs(freq, evecs)
    return freq, evecs


def logproduct_freq2(freqs, nzero, nnegative=0, eps=1e-4):
    ''' calculate the log product of positive (squared) frequencies
    
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
    '''
    inegative = 0
    izero = 0
    lnf = 0.
    n = 0
    for f in freqs:
        if np.abs(f) < eps:
            izero += 1
            continue
        if f < 0.:
            inegative += 1
            continue
        lnf += np.log(f)
        n+=1
        
    if nzero != izero:
        raise ValueError("the number of zero eigenvalues differs from the expected value")
    if nnegative != inegative:
        raise ValueError("the number of negative eigenvalues differs from the expected "
                         "number (not a minimum / transition state?)")

    return n, lnf
