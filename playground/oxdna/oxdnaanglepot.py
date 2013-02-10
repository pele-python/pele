from pygmin.utils.rbtools import CoordsAdapter
from pygmin.potentials.fortran.rmdrvt import rmdrvt
import oxdnagmin_ as GMIN
from pygmin.potentials import GMINPotential
import numpy as np

class OXDNAAnglePotential(GMINPotential):
    def __init__(self, theta0 = 140./180.*np.pi, k = 4.):
        GMINPotential.__init__(self, GMIN)
        self.theta0 = theta0
        self.k = k
        print "theta0 is", theta0 
        
    def getEnergy(self, coords):
        E = GMINPotential.getEnergy(self, coords)
        ca = CoordsAdapter(nrigid=coords.size/6, coords=coords)
        RMX = [rmdrvt(p, True) for p in ca.rotRigid]
        
        xbase = [x - 0.4*np.dot(r[0], np.array([1., 0., 0.])) for x, r in zip(ca.posRigid, RMX)]
        
        Eangle = 0
        v2 = xbase[1] - xbase[0]
        v2 /= np.linalg.norm(v2)
        for i in xrange(1,ca.nrigid-1):
            v1 = -v2.copy()
            v2 = xbase[i+1] - xbase[i]
            v2 /= np.linalg.norm(v2)
            
            theta = np.arccos(np.dot(v1,v2))
            Eangle += 0.5*self.k*(theta - self.theta0)**2
        return E + Eangle
    
    def getEnergyGradient(self, coords):
        #return self.getEnergy(coords), self.NumericalDerivative(coords)
        E, grad = GMINPotential.getEnergyGradient(self, coords)
        
        ca = CoordsAdapter(nrigid=coords.size/6, coords=coords)
        RMX = [rmdrvt(p, True) for p in ca.rotRigid]
        xbase = np.zeros([ca.nrigid, 3])
        for x, R, xb in zip(ca.posRigid, RMX, xbase):
            xb[:] = x - 0.4*np.dot(R[0], np.array([1., 0., 0.]))
        
        gbase = np.zeros_like(xbase)
        Eangle = 0.
        v2 = xbase[1] - xbase[0]
        absv2 = np.linalg.norm(v2)
        v2 /= absv2 
        for i in xrange(1,ca.nrigid-1):
            v1 = -v2
            absv1 = absv2
            v2 = xbase[i+1] - xbase[i]
            absv2 = np.linalg.norm(v2)
            v2 /= absv2
            
            v1v2 = np.dot(v1,v2)
            theta = np.arccos(v1v2)
            Eangle += 0.5*self.k*(theta - self.theta0)**2
            
            acos_prime = 1. / np.sqrt(1. - v1v2**2);
            s = self.k*(theta -  self.theta0)*acos_prime
            gbase[i-1] += s * (- v2 / absv1 +  v1v2 * v1 / absv1 )
            gbase[i]   += s * ( v1/absv2 + v2/absv1 
                              - v1v2 * v1 / absv1
                              - v1v2 * v2 / absv2)
            gbase[i+1] += s * (- v1 / absv2 +  v1v2 * v2 / absv2 )
        
        cg = CoordsAdapter(nrigid=ca.nrigid, coords=grad)    
        cg.posRigid += gbase
        for i in xrange(ca.nrigid):
            x = -0.4*np.array([1., 0., 0.])
            R = RMX[i]
            cg.rotRigid[i][0] += np.dot(gbase[i], np.dot(R[1], x))
            cg.rotRigid[i][1] += np.dot(gbase[i], np.dot(R[2], x))
            cg.rotRigid[i][2] += np.dot(gbase[i], np.dot(R[3], x))
        return E + Eangle, grad
        