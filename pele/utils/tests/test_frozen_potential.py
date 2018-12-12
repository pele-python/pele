import unittest
import numpy as np

from pele.systems import LJCluster
from pele.utils.frozen_atoms import FrozenCoordsConverter, FrozenPotWrapper, makeBLJNeighborListPotFreeze
from pele.potentials.tests import _base_test
from pele.potentials import ljpshiftfast
from pele.optimize import lbfgs_cpp


class TestFrozenCoordsConverter(unittest.TestCase):
    def setUp(self):
        self.ndof = 97
        
        self.reference_coords = np.random.rand(self.ndof)
        self.frozen_dof = np.array([0, 1, 2, 10, 40])
        
        
        self.converter = FrozenCoordsConverter(self.reference_coords, self.frozen_dof)
        
    def test1(self):
        reduced_coords = self.converter.get_reduced_coords(self.reference_coords)
        self.assertEqual(self.ndof - len(self.frozen_dof), len(reduced_coords))
        
    def test2(self):
        full_coords = self.converter.get_full_coords(self.converter.get_reduced_coords(self.reference_coords))
        self.assertTrue((full_coords == self.reference_coords).all())



_x0 = np.array([ 0.57688628,  1.63264454,  0.52232462, -0.59207205, -1.12241977,
               -0.84475259,  1.21387714,  0.22287968, -1.10619103,  1.16204747,
                1.4000165 ,  1.79543327, -1.50905004,  0.92089769,  0.07877646,
               -0.46906172,  1.67647676,  0.46896879, -0.12919939,  0.58224394,
                1.1688233 , -1.50977381, -0.29852604, -1.33057737,  1.77018274,
                1.37530175,  0.94323983, -0.63978834, -0.98402637,  0.14040699,
               -1.95924083, -0.32402232,  0.51705612, -0.04631138,  1.52063544,
               -1.4994962 ,  0.5193755 ,  0.80057462, -0.16855816,  1.5149072 ,
               -1.67154403,  0.25969931,  0.18831025, -0.56964543, -0.28784198,
                1.0201419 , -1.58242752,  1.16234538,  0.95624005, -0.58564718,
                1.43745605,  0.14717315, -1.62810577, -0.23699778, -1.25591024,
               -0.27663007, -0.30178502, -1.83432242,  0.68439946, -0.97991396,
                1.84988959, -0.35742387,  0.9532375 , -0.73335135,  0.45797021,
               -0.78473658, -1.50658925,  1.6789772 , -0.75979073,  0.59962216,
               -1.55046892, -1.16380361, -0.45807742,  0.64723181,  0.19508046,
               -0.12398379, -0.45532487,  1.29824411,  0.14173643,  1.52473984,
                1.51986249, -1.12701665,  0.52574456,  1.5239106 , -0.82181166,
               -1.9638572 , -0.20399771,  0.55839767, -0.03660041,  0.53964983,
               -1.59041491, -1.26204821, -0.19400227,  0.24893474, -0.63750162,
               -1.55903735, -0.57371507,  0.56996591, -1.77188037, -0.96037032,
               -0.20030573,  0.73928693,  0.20784111, -1.17011912,  0.66312063,
               -1.09989911, -1.56561309,  0.82926306, -0.54113531,  1.45132978,
               -0.56530919, -1.41129307,  1.24742765, -1.73101566,  1.10246955,
                1.42829452, -1.65498697,  0.46503184,  1.76533541, -0.61524605])

class TestBLJNeighborListFreeze(_base_test._TestConfiguration):
    def setUp(self):
        self.x0 = _x0
        self.e0 = -87.46393381926839
        natoms = self.x0.size // 3
        ntypeA = int(natoms*0.8)
        ntypeB = natoms - ntypeA
        freezelist = list(range(ntypeA//2)) + list(range(ntypeA,ntypeA+ntypeB//2))
        self.freezelist = freezelist
        rcut=2.5
        boxl=5.5
        self.pot = makeBLJNeighborListPotFreeze(natoms, freezelist, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
        self.blj = ljpshiftfast.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)

    def test_agains_blj(self):
        eblj = self.blj.getEnergy(self.x0)
        self.assertAlmostEqual(eblj, self.e0, places=6)
        
    def test_quench(self):
        ret = lbfgs_cpp(self.x0, self.pot)
        eblj = self.blj.getEnergy(ret.coords)
        self.assertAlmostEqual(eblj, ret.energy, places=6)
        
        # ensure the frozen atoms have not moved
        frozen_dof = self.pot.frozen1d
        self.compare_arrays(self.x0[frozen_dof], ret.coords[frozen_dof])


#class Test1(unittest.TestCase):
#    def test1(self):
#        import numpy as np
#        from pele.potentials import LJ
#        from pele.utils.frozen_atoms import FrozenPotWrapper
#        from pele.optimize import mylbfgs
#        natoms = 4
#        pot = LJ()
#        
#        reference_coords = np.random.uniform(-1, 1, [3*natoms])
#        print reference_coords
#        
#        #freeze the first two atoms (6 degrees of freedom)
#        frozen_dof = range(6)
#        
#        fpot = FrozenPotWrapper(pot, reference_coords, frozen_dof)
#        
#        reduced_coords = fpot.coords_converter.get_reduced_coords(reference_coords)
#        
#        print "the energy in the full representation:" 
#        print pot.getEnergy(reference_coords)
#        print "is the same as the energy in the reduced representation:"
#        print fpot.getEnergy(reduced_coords)
#        
#        ret = mylbfgs(reduced_coords, fpot)
#        print "after a minimization the energy is ", ret.energy, "and the rms gradient is", ret.rms
#        print "the coordinates of the frozen degrees of freedom are unchanged"
#        print "starting coords:", reference_coords
#        print "minimized coords:", fpot.coords_converter.get_full_coords(ret.coords)
#
#    def test2(self, natoms=40, boxl=5.5):
#        from pele.utils.frozen_atoms import FreezePot, makeBLJNeighborListPotFreeze
#        import pele.potentials.ljpshiftfast as ljpshift
#        from pele.optimize import mylbfgs
#        from pele.utils.neighbor_list import makeBLJNeighborListPot
#        ntypeA = int(natoms*0.8)
#        ntypeB = natoms - ntypeA
#        rcut = 2.5
#        freezelist = range(ntypeA/2) + range(ntypeA,ntypeA+ntypeB/2)
#        nfrozen = len(freezelist)
#        print "nfrozen", nfrozen
#        coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2
#        
#        
#        NLblj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)
#        blj = FreezePot(NLblj, freezelist, natoms)
#    
#        pot = makeBLJNeighborListPotFreeze(natoms, freezelist, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
#        #pot = FreezePot(NLpot, freezelist)
#
#        from pele.potentials import LJ
#        lj = LJ()
#        ret = mylbfgs(coords, lj, tol=10)        
#        ret = mylbfgs(ret.coords, pot, tol=10)
#        print repr(ret.coords)
#        print ret.energy
#        print ret.rms
#        
#        
#        eblj = blj.getEnergy(coords)
#        print "blj energy", eblj
#        
#        epot = pot.getEnergy(coords)
#        print "mcpot energy", epot
#        
#        print "difference", (epot - eblj)/eblj
#        pot.test_potential(coords)
#        print "\n"
#        
#        ret1 = mylbfgs(coords, blj, iprint=-11)
#        np.savetxt("out.coords", ret1.coords)
#        print "energy from quench1", ret1.energy
#        ret2 = mylbfgs(coords, pot, iprint=-1)
#        print "energy from quench2", ret2.energy
#        
#        print "ret1 evaluated in both potentials", pot.getEnergy(ret1.coords), blj.getEnergy(ret1.coords)
#        print "ret2 evaluated in both potentials", pot.getEnergy(ret2.coords), blj.getEnergy(ret2.coords)
#        
#        coords = ret1.coords
#        e1, g1 = blj.getEnergyGradient(coords)
#        e2, g2 = pot.getEnergyGradient(coords)
#        print "energy difference from getEnergyGradient", (e2 - e1)
#        print "largest gradient difference", np.max(np.abs(g2-g1))
#        print "rms gradients", np.linalg.norm(g1)/np.sqrt(len(g1)), np.linalg.norm(g2)/np.sqrt(len(g1))
#    
#        if True:
#            for subpot in pot.pot.potentials:
#                nl = subpot
#                print "number of times neighbor list was remade:", nl.buildcount, "out of", nl.count
#        
#        if False:
#            try: 
#                import pele.utils.pymolwrapper as pym
#                pym.start()
#                pym.draw_spheres(np.reshape(coords,[-1,3]), "A", 1)
#                pym.draw_spheres(np.reshape(ret1.coords,[-1,3]), "A", 2)
#                pym.draw_spheres(np.reshape(ret2.coords,[-1,3]), "A", 3)
#            except ImportError:
#                print "Could not draw using pymol, skipping this step" 


if __name__ == "__main__":
    unittest.main() 

