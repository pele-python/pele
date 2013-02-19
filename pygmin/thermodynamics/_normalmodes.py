import numpy as np    

__all__ = ["normalmode_frequencies", "logproduct_freq2"]

def normalmode_frequencies(hessian, metric=None):
    ''' calculate normal mode frequencies
    
    Parameters
    ----------
    hessian:
        hessian marix
        
    metric: 
        mass weighted metric tensor
    '''
    A = hessian
    if metric is not None:
        A = np.dot(np.linalg.inv(metric), hessian)
   
    return np.linalg.eigvalsh(A)

def logproduct_freq2(freqs, nzero, nnegative=0, eps=1e-4):
    ''' calculate the log product of positive frequencies
    
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
    touple of number of conisdered frequences and log product of frequencies
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