'''
Created on 13 Apr 2012

@author: ruehle
'''
import numpy as np

__all__ = ["BasePotential"]


class BasePotential(object):
    '''
    Base class for all potentials
    
    Derived classes must overload getEnergy().  It is also highly
    recommended to overload getEnergyGradient(), otherwise gradients
    will be calculated numerically  
    
        getEnergyGradient()
    '''
    def getEnergy(self, coords):
        """return the energy at the given coordinates"""
        raise NotImplementedError
        
    def getEnergyGradient(self, coords):
        """return the energy and gradient at the given coordinates"""
        return self.getEnergyGradientNumerical(coords)

    def getEnergyGradientNumerical(self, coords):
        return self.getEnergy(coords), self.NumericalDerivative(coords, 1e-8)
            
    def NumericalDerivative(self, coords, eps=1e-6):
        """return the gradient calculated numerically"""
        g = np.zeros(coords.size)
        x = coords.copy()
        for i in xrange(coords.size):
            x[i] += eps
            g[i] = self.getEnergy(x)
            x[i] -= 2. * eps
            g[i] -= self.getEnergy(x)
            g[i] /= 2. * eps
            x[i] += eps
        return g
    
    def getGradient(self, coords):
        """return the gradient at the given coordinates"""
        e, g = self.getEnergyGradient(coords)
        return g     

    #routines involving interaction lists
    def getEnergyListSlow(self, coords, ilist):
        coords2 = np.zeros(6)
        E = 0.
        for i,j  in ilist:
            coords2[:3] = coords[i*3:i*3+3]
            coords2[3:] = coords[j*3:j*3+3]
            E += self.getEnergy( coords2 )
        return E

    def getEnergyList(self, coords, ilist):
        """
        the energy of select interactions defined in ilist
        """
        return self.getEnergyListSlow(coords, ilist)

    def getEnergyGradientListSlow(self, coords, ilist):
        """you really don't want to be using this..."""
        coords2 = np.zeros(6)
        E = 0.
        g = np.zeros( len(coords) )
        for i,j  in ilist:
            coords2[:3] = coords[i*3:i*3+3]
            coords2[3:] = coords[j*3:j*3+3]
            de, dg = self.getEnergyGradient( coords2 )
            E += de
            g[i*3:i*3+3] += dg[:3]
            g[j*3:j*3+3] += dg[3:]
        return E, g

    def getEnergyGradientList(self, coords, ilist):
        """
        the energy and gradient of select interactions defined in ilist
        """
        return self.getEnergyGradientListSlow(coords, ilist)
    
    def NumericalHessian(self, coords, eps=1e-6):
        """return the Hessian matrix of second derivatives computed numerically
        
        this takes 2*len(coords) calls to getGradient
        """
        x = coords.copy()
        ndof = len(x)
        hess = np.zeros([ndof, ndof])
        for i in xrange(ndof):
            xbkup = x[i]
            x[i] += eps
            g1 = self.getGradient(x)
            x[i] = xbkup - eps
            g2 = self.getGradient(x)            
            hess[i,:] = (g1 - g2) / (2. * eps)
            x[i] = xbkup
        return hess

            
    def getEnergyGradientHessian(self, coords):
        """return the energy, gradient, and Hessian at the given coordinates"""
        e, g = self.getEnergyGradient(coords)
        hess = self.NumericalHessian(coords)
        return e, g, hess
    
    def test_potential(self, coords):
        """print some information testing whether the analytical gradients are correct"""
        E1 = self.getEnergy(coords)
        E2, grad = self.getEnergyGradient(coords)
        gradnum = self.NumericalDerivative(coords)
        print "testing energy and gradient"
        print "energy from getEnergy        ", E1
        print "energy from getEnergyGradient", E2
        print "difference", np.abs(E1-E2)
#        print "analytical gradient", grad
#        print "numerical gradient ", gradnum
        print "analytical rms gradient", np.linalg.norm(grad) / np.sqrt(coords.size)
        print "numerical rms gradient ", np.linalg.norm(gradnum) / np.sqrt(coords.size)
        print "maximum difference between analytical and numerical gradient", np.max(np.abs(grad-gradnum))
        print "normalized by the maximum gradient", np.max(np.abs(grad-gradnum)) / np.max(np.abs(grad))

class potential(BasePotential):
    """
    for backward compatibility
    """
    pass

