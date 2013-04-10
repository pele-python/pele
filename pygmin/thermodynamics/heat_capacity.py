import numpy as np


def dos_to_cv(energies, ldos, T, K=1.):
    """compute the heat capacity and other thermodynamic quantities from density of states

    Parameters
    ----------
    energies : array
        array of energies
    ldos : array
        ldos[i] is the log of the density of states at energy energies[i]
    T : array
        Temperatures at which to do the calculations.  Should be normalized by the 
        Bolzmann constant (you should pass k_B*T)
    K : int
        number of vibrational degrees of freedom (3*N-6 for clusters)
    """
    beta = np.array(1./T)

    lZ =  ldos[np.newaxis,:] - beta[:,np.newaxis] * energies[np.newaxis,:]
    
    # subtract out the smallest value to avoid overflow issues when lZ is exponentiated
    lZmax = np.max(lZ,axis=1) #  maximum lZ for each temperature
    lZ -= lZmax[:,np.newaxis]

    # compute Z, <E> and <E**2>
    Zpref = np.exp(lZ)
    Z = np.sum(Zpref, axis=1 )
    U = np.sum(Zpref * energies[np.newaxis,:], axis=1 )
    U2 = np.sum(Zpref * energies[np.newaxis,:]**2, axis=1 )
    U /= Z
    U2 /= Z

    # compute Cv from the energy fluctuations
    Cv = float(K)/2 + (U2 - U**2) * beta**2 # this is not quite right
    
    Z *= np.exp(lZmax)
    return Z, U, U2, Cv

def minima_to_cv(minima, T, k):
    """compute the heat capacity and other thermodynamic quantities from a list of minima using the harmonic approximation
    
    Parameters
    ----------
    minima : list of Minimum objects
        mimima from which to compute the thermodynamic computations
    T : array
        Temperatures at which to do the calculations.  Should be normalized by the 
        Bolzmann constant (you should pass k_B*T)
    k : int
        number of degrees of vibrational freedom (3*N-6 for clusters) 
    
    Notes
    -----
    See DJW Energy Landscapes book page 371
    Za = P * exp(-beta Ea) / (beta * h * nu_bar)**k / O_a
    """
    beta = np.array(1./T)
    k = float(k)
    Emin = min([m.energy for m in minima])
    Z = 0.
    U = 0.
    U2 = 0.
    for m in minima:
        E = m.energy
        lZpref = -beta * (E - Emin) - m.fvib/2. - np.log(m.pgorder)
        Zpref = np.exp(lZpref)
        
        Z += Zpref
        U += Zpref * E
        U2 += Zpref * E**2
    
    U /= Z
    U2 /= Z
    Cv = k + (U2 - U**2) * beta**2   # this is not completely correct
    Z *= np.exp(beta * Emin) / beta**k
    return Z, U, U2, Cv
