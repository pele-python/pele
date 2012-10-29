"""
tools for dealing with frozen atoms
"""
import numpy as np

import pygmin.potentials.ljpshift as ljpshift
from pygmin.potentials.potential import potential as basepot
from pygmin.potentials.ljcut import LJCut
from pygmin.utils.neighbor_list import NeighborListSubset, NeighborListPotential, MultiComponentSystem

__all__ = ["makeBLJNeighborListPotFreeze", "FreezePot"]


class MultiComponentSystemFreeze(basepot):
    """
    a potential wrapper for multiple potentials
    
    Parameters
    ----------
    potentials :
        a list of potential objects
    """
    def __init__(self, potentials_mobile, potentials_frozen):
        self.potentials = potentials_mobile
        self.potentials_frozen = potentials_frozen
        self.count = 0

    def _setup(self, coords):
        """
        get the energy from the frozen-frozen interactions
        """
        self.Eff = 0.
        for pot in self.potentials_frozen:
            self.Eff += pot.getEnergy(coords)

    def getEnergy(self, coords):
        if self.count == 0:
            self._setup(coords)
        self.count += 1
        E = self.Eff
        for pot in self.potentials:
            E += pot.getEnergy(coords)
        return E

    def getEnergyGradient(self, coords):
        if self.count == 0:
            self._setup(coords)
        self.count += 1
        Etot = self.Eff
        gradtot = np.zeros(np.shape(coords))
        for pot in self.potentials:
            E, grad = pot.getEnergyGradient(coords)
            Etot += E
            gradtot += grad
        return Etot, gradtot

  
def makeBLJNeighborListPotFreeze(natoms, frozenlist, ntypeA = None, rcut = 2.5, boxl=None):
    """
    create the potential object for the kob andersen binary lennard jones with frozeen particles
    
    Parameters
    ----------
    natoms :
        number of atoms in the system
    frozenlist : 
        list of frozen atoms
    ntypeA : 
        number of atoms of type A.  It is assumed that they have indices in [0,ntypeA]
    rcut : 
        the cutoff for the lj potential in units of sigA
    boxl : 
        the box length for periodic box.  None for no periodic boundary conditions
    """
    print "making BLJ neighborlist potential", natoms, ntypeA, rcut, boxl
    #rcut = 2.5
    #natoms = 40
    if ntypeA is None:
        ntypeA = int(natoms * 0.8)
    Alist = range(ntypeA)
    Blist = range(ntypeA, natoms)
    
    frozenA = np.array( [i for i in Alist if i in frozenlist] )
    mobileA = np.array( [i for i in Alist if i not in frozenlist] )
    frozenB = np.array( [i for i in Blist if i in frozenlist] )
    mobileB = np.array( [i for i in Blist if i not in frozenlist] )
    
    
    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut)

        
    ljAA = LJCut(eps=blj.AA.eps, sig=blj.AA.sig, rcut=rcut*blj.AA.sig, boxl=boxl)
    ljBB = LJCut(eps=blj.BB.eps, sig=blj.BB.sig, rcut=rcut*blj.BB.sig, boxl=boxl)
    ljAB = LJCut(eps=blj.AB.eps, sig=blj.AB.sig, rcut=rcut*blj.AB.sig, boxl=boxl)
    
    
    #ten lists in total
    #nlAA_ff
    #nlAA_mm
    #nlAA_mf
    #nlBB_ff
    #nlBB_mm
    #nlBB_mf
    #nlAB_ff
    #nlAB_mm
    #nlAB_mf
    #nlAB_fm
    
    nlAA_ff = NeighborListSubset(natoms, rcut, frozenA, boxl=boxl )
    nlAA_mm = NeighborListSubset(natoms, rcut, mobileA, boxl=boxl )
    nlAA_mf = NeighborListSubset(natoms, rcut, mobileA, Blist=frozenA, boxl=boxl )
    
    nlBB_ff = NeighborListSubset(natoms, rcut, frozenB, boxl=boxl )
    nlBB_mm = NeighborListSubset(natoms, rcut, mobileB, boxl=boxl )
    nlBB_mf = NeighborListSubset(natoms, rcut, mobileB, Blist=frozenB, boxl=boxl )
    
    nlAB_ff = NeighborListSubset(natoms, rcut, frozenA, Blist=frozenB, boxl=boxl )
    nlAB_mm = NeighborListSubset(natoms, rcut, mobileA, Blist=mobileB, boxl=boxl )
    nlAB_mf = NeighborListSubset(natoms, rcut, mobileA, Blist=frozenB, boxl=boxl )
    nlAB_fm = NeighborListSubset(natoms, rcut, mobileB, Blist=frozenA, boxl=boxl )
    
    potlist_frozen = [
                 NeighborListPotential(nlAA_ff, ljAA),
                 NeighborListPotential(nlBB_ff, ljBB),
                 NeighborListPotential(nlAB_ff, ljAB)
                ]
    potlist_mobile = [
                 NeighborListPotential(nlAA_mm, ljAA),
                 NeighborListPotential(nlAA_mf, ljAA),
                 NeighborListPotential(nlBB_mm, ljBB),
                 NeighborListPotential(nlBB_mf, ljBB),
                 NeighborListPotential(nlAB_mm, ljAB),
                 NeighborListPotential(nlAB_mf, ljAB),
                 NeighborListPotential(nlAB_fm, ljAB),
                 ]


    mcpot = MultiComponentSystemFreeze(potlist_mobile, potlist_frozen)
    frozenpot = FreezePot(mcpot, frozenlist)
    return frozenpot


class FreezePot(basepot):
    """
    potential wrapper for frozen particles
    
    the gradient will be set to zero for all frozen particles
    
    Parameters
    ----------
    pot : 
        the potential object
    frozen : 
        a list of frozen particles
    """
    def __init__(self, pot, frozen):
        self.pot = pot
        self.frozen_atoms = frozen  #a list of frozen atoms
        self.frozen1d = np.zeros(len(self.frozen_atoms)*3, np.integer) # a list of frozen coordinates (x,y,z) for each atom
        j = 0
        for i in self.frozen_atoms:
            self.frozen1d[j] = i*3
            self.frozen1d[j+1] = i*3 + 1
            self.frozen1d[j+2] = i*3 + 2
            j+=3

    def getEnergy(self, coords):
        return self.pot.getEnergy(coords)

    def getEnergyGradient(self, coords):
        e, grad = self.pot.getEnergyGradient(coords)
        grad[self.frozen1d] = 0.
        return e, grad

    def getEnergyList(self, coords):
        return self.pot.getEnergyList(coords)
        
    def getEnergyGradientList(self, coords):
        e, grad = self.pot.getEnergyGradient(coords)
        grad[self.frozen1d] = 0.
        return e, grad


#########################################################
#testing stuff below here
#########################################################

def test(natoms = 40, boxl=4.):
    import pygmin.potentials.ljpshiftfast as ljpshift
    import pygmin.defaults as defaults
    from pygmin.utils.neighbor_list import makeBLJNeighborListPot
    ntypeA = int(natoms*0.8)
    ntypeB = natoms - ntypeA
    rcut = 2.5
    freezelist = range(ntypeA/2) + range(ntypeA,ntypeA+ntypeB/2)
    nfrozen = len(freezelist)
    print "nfrozen", nfrozen
    coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2
    
    
    NLblj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)
    blj = FreezePot(NLblj, freezelist)

    pot = makeBLJNeighborListPotFreeze(natoms, freezelist, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
    #pot = FreezePot(NLpot, freezelist)
    

    
    eblj = blj.getEnergy(coords)
    print "blj energy", eblj
    
    epot = pot.getEnergy(coords)
    print "mcpot energy", epot
    
    print "difference", (epot - eblj)/eblj
    
    ret1 = defaults.quenchRoutine(coords, blj.getEnergyGradient, iprint=-11)
    np.savetxt("out.coords", ret1[0])
    print "energy from quench1", ret1[1]
    ret2 = defaults.quenchRoutine(coords, pot.getEnergyGradient, iprint=-1)
    print "energy from quench2", ret2[1]
    
    print "ret1 evaluated in both potentials", pot.getEnergy(ret1[0]), blj.getEnergy(ret1[0])
    print "ret2 evaluated in both potentials", pot.getEnergy(ret2[0]), blj.getEnergy(ret2[0])
    
    coords = ret1[0]
    e1, g1 = blj.getEnergyGradient(coords)
    e2, g2 = pot.getEnergyGradient(coords)
    print "energy difference from getEnergyGradient", (e2 - e1)
    print "largest gradient difference", np.max(np.abs(g2-g1))
    print "rms gradients", np.linalg.norm(g1)/np.sqrt(len(g1)), np.linalg.norm(g2)/np.sqrt(len(g1))

    
    for subpot in pot.pot.potentials:
        nl = subpot.neighborList
        print "number of times neighbor list was remade", nl.buildcount, "out of", nl.count
    
    if False:
        try: 
            import pygmin.utils.pymolwrapper as pym
            pym.start()
            pym.draw_spheres(np.reshape(coords,[-1,3]), "A", 1)
            pym.draw_spheres(np.reshape(ret1[0],[-1,3]), "A", 2)
            pym.draw_spheres(np.reshape(ret2[0],[-1,3]), "A", 3)
        except ImportError:
            print "Could not draw using pymol, skipping this step" 

if __name__ == "__main__":
    test()