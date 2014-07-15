from pele.utils.rbtools import CoordsAdapter
from pele.angleaxis.aatopology import rotMatDeriv
import oxdnagmin_ as GMIN
from pele.potentials import GMINPotential
import numpy as np
from dihedral import dihedral_angle, dihedral_gradient

def U_torsion_back(theta):
    return 0.5*theta**2

def U_torsion_back_grad(theta):
    return 0.5*theta**2, theta

class OXDNAAnglePotential(GMINPotential):
    def __init__(self, theta0 = 140./180.*np.pi, k = 4., use_torsion=False):
        GMINPotential.__init__(self, GMIN)
        self.theta0 = theta0
        self.k = k
        print "theta0 is", theta0
        self.use_torsion = use_torsion 
        
    def getEnergy(self, coords):
        E = GMINPotential.getEnergy(self, coords)
        ca = CoordsAdapter(nrigid=coords.size/6, coords=coords)
        RMX = [rotMatDeriv(p, True) for p in ca.rotRigid]
        
        xback = np.array([x - 0.4*np.dot(r[0], np.array([1., 0., 0.])) for x, r in zip(ca.posRigid, RMX)])
        xbase = np.array([x + 0.4*np.dot(r[0], np.array([1., 0., 0.])) for x, r in zip(ca.posRigid, RMX)])
        
        Eangle = 0
        v2 = xback[1] - xback[0]
        v2 /= np.linalg.norm(v2)
        for i in xrange(1,ca.nrigid-1):
            v1 = -v2.copy()
            v2 = xback[i+1] - xback[i]
            v2 /= np.linalg.norm(v2)
            
            theta = np.arccos(np.dot(v1,v2))
            Eangle += 0.5*self.k*(theta - self.theta0)**2
            
        # add the torsion angle
        Etorsion = 0
        if self.use_torsion:
            for i in xrange(ca.nrigid-3):
                theta = dihedral_angle(xback[i:i+4])
                Etorsion += U_torsion_back(theta)
        return E + Eangle + Etorsion
    
    def getEnergyGradient(self, coords):
        #return self.getEnergy(coords), self.NumericalDerivative(coords)
        E, grad = GMINPotential.getEnergyGradient(self, coords)
        
        ca = CoordsAdapter(nrigid=coords.size/6, coords=coords)
        RMX = [rotMatDeriv(p, True) for p in ca.rotRigid]

        xback = [x - 0.4*np.dot(r[0], np.array([1., 0., 0.])) for x, r in zip(ca.posRigid, RMX)]
        xbase = [x - 0.4*np.dot(r[0], np.array([1., 0., 0.])) for x, r in zip(ca.posRigid, RMX)]
        
        gback = np.zeros_like(xback)
        Eangle = 0.
        v2 = xback[1] - xback[0]
        absv2 = np.linalg.norm(v2)
        v2 /= absv2 
        for i in xrange(1,ca.nrigid-1):
            v1 = -v2
            absv1 = absv2
            v2 = xback[i+1] - xback[i]
            absv2 = np.linalg.norm(v2)
            v2 /= absv2
            
            v1v2 = np.dot(v1,v2)
            theta = np.arccos(v1v2)
            Eangle += 0.5*self.k*(theta - self.theta0)**2
            
            acos_prime = 1. / np.sqrt(1. - v1v2**2);
            s = self.k*(theta -  self.theta0)*acos_prime
            gback[i-1] += s * (- v2 / absv1 +  v1v2 * v1 / absv1 )
            gback[i]   += s * ( v1/absv2 + v2/absv1 
                              - v1v2 * v1 / absv1
                              - v1v2 * v2 / absv2)
            gback[i+1] += s * (- v1 / absv2 +  v1v2 * v2 / absv2 )
        
        Etorsion = 0
        if self.use_torsion:
            for i in xrange(ca.nrigid-3):
                r = xback[i:i+4]
                theta = dihedral_angle(r)
                e_theta, g_theta = U_torsion_back_grad(theta)
                Etorsion += e_theta        
                gback[i:i+4] += g_theta * dihedral_gradient(r)  

        
        cg = CoordsAdapter(nrigid=ca.nrigid, coords=grad)    
        cg.posRigid += gback
        for i in xrange(ca.nrigid):
            x = -0.4*np.array([1., 0., 0.])
            R = RMX[i]
            cg.rotRigid[i][0] += np.dot(gback[i], np.dot(R[1], x))
            cg.rotRigid[i][1] += np.dot(gback[i], np.dot(R[2], x))
            cg.rotRigid[i][2] += np.dot(gback[i], np.dot(R[3], x))
        return E + Eangle + Etorsion, grad
        
