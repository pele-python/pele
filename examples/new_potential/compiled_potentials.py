import numpy as np
from mypotential import MyPot


def compare_potentials(pot1, pot2, coords):
    e1, g1 = pot1.getEnergyGradient(coords)
    e2, g2 = pot2.getEnergyGradient(coords)
    
    ediff = (e1-e2)
    assert ediff < 1e-5
    
    max_gdiff = np.max(np.abs(g1 - g2))
    assert max_gdiff < 1e-5

if __name__ == "__main__":
    
    natoms = 13
    pot1 = MyPot(natoms)
    e = 101.
    while e > 100:
        coords = np.random.uniform(-5,5,[natoms*3])
        e = pot1.getEnergy(coords)
    
    try:
        from fortran_potential.mypotential import MyPotFortran
        pot2 = MyPotFortran(natoms)
        print "testing fortran potential"
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import fortran potential"