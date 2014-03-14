import numpy as np
from mypotential import MyPot


def compare_potentials(pot1, pot2, coords):
    e1, g1 = pot1.getEnergyGradient(coords)
    e2, g2 = pot2.getEnergyGradient(coords)
    
    ediff = (e1-e2)
    assert ediff < 1e-5
    
    max_gdiff = np.max(np.abs(g1 - g2))
    if max_gdiff > 1e-5:
        print "maximum difference in the gradient", max_gdiff
        print e1, e2
        print g1
        print g2

if __name__ == "__main__":
    
    natoms = 13
    pot1 = MyPot(natoms)
    e = 1001.
    while e > 1000:
        coords = np.random.uniform(-3,3,[natoms*3])
        e = pot1.getEnergy(coords)
    
    try:
        from fortran_potential.mypotential import MyPotFortran
        pot2 = MyPotFortran(natoms)
        print "testing fortran potential"
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import fortran potential"

    try:
        from cython_potential.mypotential import MyPotCython
        pot2 = MyPotCython(natoms)
        print "testing cython potential"
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import cython potential"

    try:
        from cpp_potential.mypotential import MyPotCpp
        pot2 = MyPotCpp(natoms)
        print "testing c++ potential"
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import c++ potential"
    try:
        from c_potential.mypotential import MyPotC
        pot2 = MyPotC(natoms)
        print "testing c potential"
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import c potential"
    try:
        from fortran_iso_c_bindings.mypotential import MyPotIsoC
        print "testing fortran potential with iso_c_bindings"
        pot2 = MyPotIsoC(natoms)
        compare_potentials(pot1, pot2, coords)
    except ImportError:
        print "could not import fortran iso_c_bindings potential"
