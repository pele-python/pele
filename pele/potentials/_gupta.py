import numpy as np

from pele.potentials import BasePotential
from pele.potentials.fortran import gupta as fortran_gupta

class gupta(BasePotential):
    """The Gupta potential for transition metals.

    Parameters
    ----------
     p : float
         p is dimensionless, goes inside the repulsive exponential
     q : float
         q is dimensionless, goes inside the attractive exponential 
         (i.e. the pseudo-density)
     A : float
         A is an energy pre-factor for the pairwise repulsion
    xi : float
         xi is an energy pre-factor for the many-body attraction

    Notes
    -----
    The default parameter set is for Platinum, with the values 
    fitted to bulk properties and also adequate for clusters.

    Here is a selection of representative parameter values:
    [ Cleri and Rosato, PRB 48, 22 (1993) ]

    Metal      p       q     A (eV)  xi (eV)
    
      Ni    16.999   1.189   0.0376   1.070
      Cu    10.960   2.278   0.0855   1.224
      Rh    18.450   1.867   0.0629   1.660  
      Pd    10.867   3.742   0.1746   1.718
      Ag    10.928   3.139   0.1028   1.178
      Ir    16.980   2.691   0.1156   2.289
      Pt    10.612   4.004   0.2975   2.695
      Au    10.229   4.036   0.2061   1.790
      Al     8.612   2.516   0.1221   1.316
      Pb     9.576   3.648   0.0980   0.914
      V      5.206   1.220   0.6124   2.441
    
      V*     6.82    1.85    0.10     1.88   [Phys. B 215, 377 (1995)]
      V^     6.3579  1.6709  0.2469   2.0396 [PSS B 215, 1127 (1999)]

    """
    def __init__(self, p=10.612, q=4.004, A=0.2975, xi=2.695):
        self.p = p
        self.q = q
        self.A = A
        self.xi = xi
        
    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e

    def getEnergyGradient(self, coords):
        g, e = fortran_gupta.gupta(coords, self.p, self.q, 
                                   self.A, self.xi)
        return e, g

#
# testing only below here
#
def test_platinum():
    pot = gupta()
    natoms = 20
    x = np.random.uniform(-2,2,[3*natoms])
    pot.getEnergy(x)
    pot.test_potential(x)

if __name__ == "__main__":
    test_platinum()

