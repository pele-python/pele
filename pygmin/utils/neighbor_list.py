import numpy as np
from pygmin.potentials.potential import potential as basepot
import _fortran_utils

__all__ = ["NeighborList", "NeighborListSubset", "NeighborListPotential", "MultiComponentSystem", 
           "makeBLJNeighborListPot"]

class NeighborList(object):
    """
    this class will create and keep a neighbor list updated
    """
    buildcount = 0
    def __init__(self, natoms, rcut, rskin = 0.5, boxl = None):
        """
        :rcut: the cutoff distance for the potential
        
        :rskin: the skin distance.  atoms are listed as 
            neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
            will have more interaction pairs, but will need to be updated less
            frequently
        """
        self.oldcoords = np.zeros([natoms,3])
        self.rcut = rcut
        self.rskin = rskin
        self.redo_displacement = self.rskin / 2.
        self.rlist = self.rcut + self.rskin
        self.rlist2 = self.rlist**2
        
        #natoms = len(coords) / 3
        self.neib_list = np.array(natoms*(natoms-1), 2)
        
        #self.buildList(coords)
    
    def buildList(self, coords):
        """
        return a list of neighbor pairs
        """
        self.buildcount += 1
        coords = np.reshape(coords, [-1,3])
        self.oldcoords = np.copy(coords)
        natoms = len(coords)/3
        nlist = 0
        for i in range(natoms):
            for j in range(i):
                R2 = sum((coords[i,:] - coords[j,:])**2)
                if R2 <= self.rlist2:
                    self.neib_list[nlist,0] = i 
                    self.neib_list[nlist,1] = j
                    nlist += 1
        self.nlist = nlist
    
    def needNewList(self, coords):
        coords = np.reshape(coords, [-1,3])
        maxR2 = np.max( ((coords - self.oldcoords)**2).sum(1) )
        return maxR2 > self.redo_displacement**2

    def getList(self, coords):
        if self.needNewList(coords):
            self.buildList(coords)
        return self.neib_list[:self.nlist]
                
            
class NeighborListSubset(object):
    """
    this class will create and keep a neighbor list updated
    """
    buildcount = 0
    count = 0
    def __init__(self, natoms, rcut, Alist, Blist = None, rskin = 0.5, boxl = None):
        """
        :rcut: the cutoff distance for the potential
        
        :rskin: the skin distance.  atoms are listed as 
            neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
            will have more interaction pairs, but will need to be updated less
            frequently
        
        :Alist: the list of atoms that are interacting
        
        :Blist: the list of atoms that are interacting with Alist
            if Blist is None or Blist is Alist then the atoms in Alist will be 
            assumed to be interacting with each other.  Duplicate interactions will
            be avoided.
        """
        self.rcut = rcut
        self.rskin = rskin
        self.redo_displacement = self.rskin / 2.
        self.rlist = self.rcut + self.rskin
        self.rlist2 = self.rlist**2
        
        if boxl is None:
            self.periodic = False
        else:
            self.periodic = True
        self.boxl = boxl
        
        self.Alist = np.array(np.copy(Alist))
        if Blist is None or Blist is Alist:
            self.onelist = True
            self.Blist = None
        else:
            self.onelist = False
            self.Blist = np.array(np.copy(Blist))
        #print "onelist", self.onelist

        if self.onelist:
            listmaxlen = len(self.Alist)*(len(self.Alist)-1)/2
        else:
            listmaxlen = len(self.Alist)*len(self.Blist)
        self.neib_list = np.zeros([listmaxlen, 2], np.integer)
        self.nlistmax = listmaxlen
        self.nlist = 0
        #print "shape neib_list", np.shape(self.neib_list)
        
        self.oldcoords = np.zeros([natoms,3])            
        #self.buildList(coords)
    
    def buildList(self, coords):
        #neib_list = np.reshape(self.neib_list, -1)
        self.buildcount += 1
        self.oldcoords = np.copy(np.reshape(coords,[-1,3]))
        if self.onelist:
            #nlist = _fortran_utils.build_neighbor_list1(
            #        coords, self.Alist, neib_list, self.rlist2)
            if self.periodic:
                neib_list, nlist = _fortran_utils.build_neighbor_list1_periodic(
                        coords, self.Alist, self.nlistmax*2, self.rlist2, self.boxl)
            else:
                neib_list, nlist = _fortran_utils.build_neighbor_list1(
                        coords, self.Alist, self.nlistmax*2, self.rlist2)
        else:
            if self.periodic:
                neib_list, nlist = _fortran_utils.build_neighbor_list2_periodic(
                        coords, self.Alist, self.Blist, self.nlistmax*2, 
                        self.rlist2, self.boxl)
            else:
                neib_list, nlist = _fortran_utils.build_neighbor_list2(
                        coords, self.Alist, self.Blist, self.nlistmax*2, self.rlist2)
        self.neib_list = np.reshape(neib_list, [-1,2])
        self.nlist = nlist
        #self.neib_list[:self.nlist,:] -= 1 #convert from fortran indices
        #print "nlist from fortran", nlist, self.neib_list[0,:], self.neib_list[self.nlist-1,:]
      
    
    def buildListSlow(self, coords):
        """
        return a list of neighbor pairs
        """
        self.buildcount += 1
        if not self.periodic:
            self.buildListFortran(coords)
            #return
        #print "rebuilding Neighbor List", self.buildcount
        coords = np.reshape(coords, [-1,3])
        self.oldcoords = np.copy(coords)
        nlist = 0
        if self.onelist:
            for k1 in range(len(self.Alist)):
                i = self.Alist[k1]
                for k2 in range(k1):
                    j = self.Alist[k2]
                    if self.periodic:
                        dr = (coords[i,:] - coords[j,:])
                        dr -= self.boxl * np.round( dr / self.boxl )
                        R2 = np.sum( dr**2 )
                    else:
                        R2 = np.sum((coords[i,:] - coords[j,:])**2)
                    if R2 <= self.rlist2:
                        #print nlist, np.shape(self.neib_list)
                        self.neib_list[nlist,0] = i 
                        self.neib_list[nlist,1] = j
                        nlist += 1
        else:
            for i in self.Alist:
                for j in self.Blist:
                    if self.periodic:
                        dr = (coords[i,:] - coords[j,:])
                        dr -= self.boxl * np.round( dr / self.boxl )
                        R2 = np.sum( dr**2 )
                    else:
                        R2 = np.sum((coords[i,:] - coords[j,:])**2)
                    if R2 <= self.rlist2:
                        self.neib_list[nlist,0] = i 
                        self.neib_list[nlist,1] = j
                        nlist += 1
        self.nlist = nlist
        #print "rebuild list: nlist", nlist
        #if not self.periodic:
        #    print "nlist not from fortran", nlist, self.neib_list[0,:], self.neib_list[self.nlist-1,:]

    
    def needNewList(self, coords):
        """
        check if any atom has moved far enough that we need to redo the neighbor list
        """
        coords = np.reshape(coords, [-1,3])
        maxR2 = np.max( ((coords[self.Alist,:] 
                          - self.oldcoords[self.Alist,:])**2).sum(1) )
        if not self.onelist:    
                tempmaxR2 = np.max( ((coords[self.Blist,:] 
                                     - self.oldcoords[self.Blist,:])**2).sum(1) )
                maxR2 = max([maxR2, tempmaxR2])
        #print "maxR2", maxR2, self.redo_displacement**2
        return maxR2 > self.redo_displacement**2

    def getList(self, coords):
        self.count += 1
        if self.needNewList(coords):
            self.buildList(coords)
        return self.neib_list[:self.nlist,:]

class NeighborListPotential(basepot):
    """
    a potential wrapper for a neighbor list
    """
    def __init__(self, neighborList, pot):
        self.neighborList = neighborList
        self.pot = pot

    def getEnergy(self, coords):
        list = self.neighborList.getList(coords)
        return self.pot.getEnergyList(coords, list)
    def getEnergyGradient(self, coords):
        list = self.neighborList.getList(coords)
        return self.pot.getEnergyGradientList(coords, list)


class MultiComponentSystem(basepot):
    """
    a potential wrapper for multiple potentials
    """
    def __init__(self, potentials):
        self.potentials = potentials
    def getEnergy(self, coords):
        E = 0.
        for pot in self.potentials:
            E += pot.getEnergy(coords)
        return E
    def getEnergyGradient(self, coords):
        Etot = 0.
        gradtot = np.zeros(np.shape(coords))
        for pot in self.potentials:
            E, grad = pot.getEnergyGradient(coords)
            Etot += E
            gradtot += grad
        return Etot, gradtot


def makeBLJNeighborListPot(natoms, ntypeA = None, rcut = 2.5, boxl=None):
    """
    recreate the binary lj with atom typea A,B from 3 interaction lists AA, BB, AB
    """
    from pygmin.potentials.ljcut import LJCut as LJ
    import pygmin.potentials.ljpshiftfast as ljpshift
    #rcut = 2.5
    #natoms = 40
    if ntypeA is None:
        ntypeA = int(natoms * 0.8)
    Alist = range(ntypeA)
    Blist = range(ntypeA, natoms)
    
    
    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut)

        
    ljAA = LJ(eps=blj.AA.eps, sig=blj.AA.sig, rcut=rcut*blj.AA.sig, boxl=boxl)
    ljBB = LJ(eps=blj.BB.eps, sig=blj.BB.sig, rcut=rcut*blj.AA.sig, boxl=boxl)
    ljAB = LJ(eps=blj.AB.eps, sig=blj.AB.sig, rcut=rcut*blj.AA.sig, boxl=boxl)
    
    nlAA = NeighborListSubset(natoms, rcut, Alist, boxl=boxl )
    nlBB = NeighborListSubset(natoms, rcut, Blist, boxl=boxl )
    nlAB = NeighborListSubset(natoms, rcut, Alist, Blist, boxl=boxl)
    
    potlist = [ 
               NeighborListPotential(nlAA, ljAA),
               NeighborListPotential(nlBB, ljBB),
               NeighborListPotential(nlAB, ljAB)
                ]
    mcpot = MultiComponentSystem(potlist)
    return mcpot

def test(natoms = 40, boxl=None):
    import pygmin.potentials.ljpshiftfast as ljpshift
    import pygmin.defaults as defaults
    ntypeA = int(natoms*0.8)
    rcut = 2.5
    coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2
    
    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)

    pot = makeBLJNeighborListPot(natoms, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
    

    
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

    
    for subpot in pot.potentials:
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
    
def test2():
    import pygmin.potentials.ljpshiftfast as ljpshift
    import pygmin.defaults as defaults
    fname = "/scratch/scratch2/js850/library/cluster/spherical/1620/PTMC/q4/oneatom/cavity200-8/ts/coords1.quench"
    fname = "out.coords"
    if False:
        coords = np.array(np.loadtxt(fname))
    else:
        natoms = 200
        coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2

    natoms = len(coords) /3
    ntypeA = int(natoms*0.8)
    rcut = 2.5
    boxl = 4
    
    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)

    pot = makeBLJNeighborListPot(natoms, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
    
    eblj = blj.getEnergy(coords)
    print "blj energy", eblj
    
    epot = pot.getEnergy(coords)
    print "mcpot energy", epot
    
    print "energy difference", (epot - eblj)
    
    e1, g1 = blj.getEnergyGradient(coords)
    e2, g2 = pot.getEnergyGradient(coords)
    print "energy difference from getEnergyGradient", (e2 - e1)
    print "largest gradient difference", np.max(np.abs(g2-g1))
    print "rms gradients", np.linalg.norm(g1)/np.sqrt(len(g1)), np.linalg.norm(g2)/np.sqrt(len(g1))


    
if __name__ == "__main__":
    #test2()
    test(natoms=100, boxl=None)
    
       
        
