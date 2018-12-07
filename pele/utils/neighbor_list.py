"""
classes to build and maintain neighborlists
.. currentmodule:: pele.utils.neighbor_list

.. autosummary:: 
    :toctree: generated/
    
    MultiComponentSystem
    NeighborListSubsetBuild
    NeighborListPotentialBuild
    NeighborListPotentialMulti

    
"""

import numpy as np

from pele.potentials.potential import potential as basepot
from . import _fortran_utils
from pele.potentials.ljcut import LJCut as LJ
import pele.potentials.ljpshiftfast as ljpshift


__all__ = ["MultiComponentSystem", "NeighborListSubsetBuild", "NeighborListPotentialBuild",
           "NeighborListPotentialMulti"]

# class NeighborList(object):
# """
# Create a neighbor list and keep it updated
#    
#    Parameters
#    ----------
#    rcut : 
#        the cutoff distance for the potential
#    rskin : 
#        the skin distance.  atoms are listed as 
#        neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
#        will have more interaction pairs, but will need to be updated less
#        frequently
#    boxl : 
#        if not None, then the system is in a periodic box of size boxl
#
#    """
#    def __init__(self, natoms, rcut, rskin = 0.5, boxl = None):
#        self.buildcount = 0
#        self.oldcoords = np.zeros([natoms,3])
#        self.rcut = rcut
#        self.rskin = rskin
#        self.redo_displacement = self.rskin / 2.
#        self.rlist = self.rcut + self.rskin
#        self.rlist2 = self.rlist**2
#        
#        #natoms = len(coords) / 3
#        self.neib_list = np.array(natoms*(natoms-1), 2)
#        
#        #self.buildList(coords)
#    
#    def buildList(self, coords):
#        """
#        return a list of neighbor pairs
#        """
#        self.buildcount += 1
#        coords = np.reshape(coords, [-1,3])
#        self.oldcoords = np.copy(coords)
#        natoms = len(coords)/3
#        nlist = 0
#        for i in range(natoms):
#            for j in range(i):
#                R2 = sum((coords[i,:] - coords[j,:])**2)
#                if R2 <= self.rlist2:
#                    self.neib_list[nlist,0] = i 
#                    self.neib_list[nlist,1] = j
#                    nlist += 1
#        self.nlist = nlist
#    
#    def needNewList(self, coords):
#        coords = np.reshape(coords, [-1,3])
#        maxR2 = np.max( ((coords - self.oldcoords)**2).sum(1) )
#        return maxR2 > self.redo_displacement**2
#
#    def getList(self, coords):
#        if self.needNewList(coords):
#            self.buildList(coords)
#        return self.neib_list[:self.nlist]


#class NeighborListSubset(object):
#    """
#    Create a neighbor list and keep it updated for a subset of the atoms.
#    
#    This class is designed to deal with only a subset of all atoms
#
#    Parameters
#    ----------
#    natoms :
#        number of atoms
#    rcut : 
#        the cutoff distance for the potential
#    rskin : 
#        the skin distance.  atoms are listed as 
#        neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
#        will have more interaction pairs, but will need to be updated less
#        frequently
#    Alist : 
#        The list of atoms that are interacting
#    Blist : the list of atoms that are interacting with Alist
#        if Blist is None or Blist is Alist then the atoms in Alist will be 
#        assumed to be interacting with each other.  Duplicate interactions will
#        be avoided.
#    boxl : 
#        if not None, then the system is in a periodic box of size boxl
#    """
#    def __init__(self, natoms, rcut, Alist, Blist = None, rskin = 0.5, boxl = None):
#        self.buildcount = 0
#        self.count = 0
#        self.rcut = rcut
#        self.rskin = rskin
#        self.redo_displacement = self.rskin / 2.
#        self.rlist = self.rcut + self.rskin
#        self.rlist2 = self.rlist**2
#        
#        if boxl is None:
#            self.periodic = False
#        else:
#            self.periodic = True
#        self.boxl = boxl
#        
#        self.Alist = np.array(np.copy(Alist), np.int64)
#        if Blist is None or Blist is Alist:
#            self.onelist = True
#            self.Blist = None
#        else:
#            self.onelist = False
#            self.Blist = np.array(np.copy(Blist))
#        #print "onelist", self.onelist
#
#        if self.onelist:
#            listmaxlen = len(self.Alist)*(len(self.Alist)-1)/2
#        else:
#            listmaxlen = len(self.Alist)*len(self.Blist)
#        self.neib_list = np.zeros([listmaxlen, 2], np.integer)
#        self.nlistmax = listmaxlen
#        self.nlist = 0
#        #print "shape neib_list", np.shape(self.neib_list)
#        
#        self.oldcoords = np.zeros([natoms,3])            
#        #self.buildList(coords)
#        
#        if self.onelist:
#            self.atomlist = list(self.Alist)
#        else:
#            self.atomlist = list(self.Alist) + list(self.Blist)
#        self.atomlist = sorted(self.atomlist)
#        
#        #we must specify the type of integer so that we can
#        #pass it to fortran without copying
#        self.atomlist = np.array(self.atomlist, np.int64)
#        self.Alist = np.array(self.Alist, np.int64)
#        if not self.onelist:
#            self.Blist = np.array(self.Blist, np.int64)
#    
#    def buildList(self, coords):
#        #neib_list = np.reshape(self.neib_list, -1)
#        self.buildcount += 1
#        self.oldcoords = np.copy(np.reshape(coords,[-1,3]))
##        raw_input("press enter to continue: onelist %d, len(alist)=%d, len(coords)=%d" % (self.onelist, len(self.Alist), len(coords)))
#        if self.onelist:
#            #nlist = _fortran_utils.build_neighbor_list1(
#            #        coords, self.Alist, neib_list, self.rlist2)
#            if self.periodic:
#                neib_list, nlist = _fortran_utils.build_neighbor_list1_periodic(
#                        coords, self.Alist, self.nlistmax*2, self.rlist2, self.boxl)
#            else:
#                neib_list, nlist = _fortran_utils.build_neighbor_list1(
#                        coords, self.Alist, self.nlistmax*2, self.rlist2)
#        else:
#            if self.periodic:
#                neib_list, nlist = _fortran_utils.build_neighbor_list2_periodic(
#                        coords, self.Alist, self.Blist, self.nlistmax*2, 
#                        self.rlist2, self.boxl)
#            else:
#                neib_list, nlist = _fortran_utils.build_neighbor_list2(
#                        coords, self.Alist, self.Blist, self.nlistmax*2, self.rlist2)
#        self.neib_list = np.reshape(neib_list, [-1,2])
#        self.nlist = nlist
#      
#    
#    def buildListSlow(self, coords):
#        """
#        return a list of neighbor pairs
#        """
#        self.buildcount += 1
#        if not self.periodic:
#            self.buildListFortran(coords)
#            #return
#        #print "rebuilding Neighbor List", self.buildcount
#        coords = np.reshape(coords, [-1,3])
#        self.oldcoords = np.copy(coords)
#        nlist = 0
#        if self.onelist:
#            for k1 in range(len(self.Alist)):
#                i = self.Alist[k1]
#                for k2 in range(k1):
#                    j = self.Alist[k2]
#                    if self.periodic:
#                        dr = (coords[i,:] - coords[j,:])
#                        dr -= self.boxl * np.round( dr / self.boxl )
#                        R2 = np.sum( dr**2 )
#                    else:
#                        R2 = np.sum((coords[i,:] - coords[j,:])**2)
#                    if R2 <= self.rlist2:
#                        #print nlist, np.shape(self.neib_list)
#                        self.neib_list[nlist,0] = i 
#                        self.neib_list[nlist,1] = j
#                        nlist += 1
#        else:
#            for i in self.Alist:
#                for j in self.Blist:
#                    if self.periodic:
#                        dr = (coords[i,:] - coords[j,:])
#                        dr -= self.boxl * np.round( dr / self.boxl )
#                        R2 = np.sum( dr**2 )
#                    else:
#                        R2 = np.sum((coords[i,:] - coords[j,:])**2)
#                    if R2 <= self.rlist2:
#                        self.neib_list[nlist,0] = i 
#                        self.neib_list[nlist,1] = j
#                        nlist += 1
#        self.nlist = nlist
#        #print "rebuild list: nlist", nlist
#        #if not self.periodic:
#        #    print "nlist not from fortran", nlist, self.neib_list[0,:], self.neib_list[self.nlist-1,:]
#
#    
#    def needNewList(self, coords):
#        """
#        check if any atom has moved far enough that we need to redo the neighbor list
#        """
#        oldcoords = self.oldcoords.reshape(-1)
#        boxl = self.boxl
#        if boxl is None:
#            boxl = 1.
##        raw_input("press enter to continue: onelist %d, len(atomlist)=%d, len(coords)=%d" % (self.onelist, len(self.atomlist), len(coords)))
#        rebuild = _fortran_utils.check_neighbor_lists(oldcoords, coords, self.atomlist,
#                                                      self.redo_displacement, self.periodic, 
#                                                      boxl)
#        rebuild = bool(rebuild)
#        if False:
#            #testing
#            rebuild_alt = self.needNewListSlow(coords)
#            if rebuild != rebuild_alt:
#                print "rebuild is incorrect", rebuild, rebuild_alt, self.periodic, self.onelist
#                print "    ", self.atomlist[-3:]
#                
#        return rebuild
#    
#    def needNewListSlow(self, coords):
#        """
#        check if any atom has moved far enough that we need to redo the neighbor list
#        """
#        coords = np.reshape(coords, [-1,3])
#        maxR2 = np.max( ((coords[self.Alist,:] 
#                          - self.oldcoords[self.Alist,:])**2).sum(1) )
#        if not self.onelist:    
#                tempmaxR2 = np.max( ((coords[self.Blist,:] 
#                                     - self.oldcoords[self.Blist,:])**2).sum(1) )
#                maxR2 = max([maxR2, tempmaxR2])
#        #print "maxR2", maxR2, self.redo_displacement**2
#        return maxR2 > self.redo_displacement**2
#
#    def getList(self, coords):
#        self.count += 1
#        if self.needNewList(coords):
#            self.buildList(coords)
#        return self.neib_list[:self.nlist,:]


#class NeighborListPotential(basepot):
#    """
#    a potential wrapper for a neighbor list
#    
#    Parameters
#    ----------
#    neighborList : 
#        the neighbor list object
#    pot :
#        the potential object
#    """
#    def __init__(self, neighborList, pot):
#        self.neighborList = neighborList
#        self.pot = pot
#
#    def getEnergy(self, coords):
#        list = self.neighborList.getList(coords)
#        return self.pot.getEnergyList(coords, list)
#    def getEnergyGradient(self, coords):
#        list = self.neighborList.getList(coords)
#        return self.pot.getEnergyGradientList(coords, list)


class NeighborListSubsetBuild(basepot):
    """
    The same as NeighborListSubset except only do the building.
    
    This class will build the neighbor lists, but can't
    check whether they need to be rebuilt.  This is meant to be used
    in situations where multiple neighbor lists are being maintained and
    we don't want to do repetative checks for whether the lists should be rebuilt
    
    It has the added benefit that no copy of coords is saved, so it takes up less memory. 

    Parameters
    ----------
    natoms :
        number of atoms
    rcut : 
        the cutoff distance for the potential
    rskin : 
        the skin distance.  atoms are listed as 
        neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
        will have more interaction pairs, but will need to be updated less
        frequently
    Alist : 
        The list of atoms that are interacting
    Blist : the list of atoms that are interacting with Alist
        if Blist is None or Blist is Alist then the atoms in Alist will be 
        assumed to be interacting with each other.  Duplicate interactions will
        be avoided.
    boxl : 
        if not None, then the system is in a periodic box of size boxl
    """

    def __init__(self, natoms, rcut, Alist, Blist=None, rskin=0.5, boxl=None):
        self.natoms = natoms
        self.buildcount = 0
        self.count = 0
        self.rcut = rcut
        self.rskin = rskin
        self.redo_displacement = self.rskin / 2.
        self.rlist = self.rcut + self.rskin
        self.rlist2 = self.rlist ** 2

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

        if self.onelist:
            listmaxlen = len(self.Alist) * (len(self.Alist) - 1) // 2
        else:
            listmaxlen = len(self.Alist) * len(self.Blist)
        #self.neib_list = np.zeros([listmaxlen, 2], np.integer)
        self.nlistmax = listmaxlen
        #self.nlist = 0
        #print "shape neib_list", np.shape(self.neib_list)

        #we must specify the type of integer so that we can
        #pass it to fortran without copying
        self.Alist = np.array(self.Alist, np.int64)
        if not self.onelist:
            self.Blist = np.array(self.Blist, np.int64)


    def buildList(self, coords):
        #neib_list = np.reshape(self.neib_list, -1)
        self.buildcount += 1
        if self.onelist:
            #nlist = _fortran_utils.build_neighbor_list1(
            #        coords, self.Alist, neib_list, self.rlist2)
            if self.periodic:
                neib_list, nlist = _fortran_utils.build_neighbor_list1_periodic(
                    coords, self.Alist, self.nlistmax * 2, self.rlist2, self.boxl)
            else:
                neib_list, nlist = _fortran_utils.build_neighbor_list1(
                    coords, self.Alist, self.nlistmax * 2, self.rlist2)
        else:
            if self.periodic:
                neib_list, nlist = _fortran_utils.build_neighbor_list2_periodic(
                    coords, self.Alist, self.Blist, self.nlistmax * 2,
                    self.rlist2, self.boxl)
            else:
                neib_list, nlist = _fortran_utils.build_neighbor_list2(
                    coords, self.Alist, self.Blist, self.nlistmax * 2, self.rlist2)
        neib_list = np.reshape(neib_list, [-1, 2])
        return neib_list[:nlist, :]


class NeighborListPotentialBuild(basepot):
    """
    a potential wrapper for a neighbor list, but only rebuild when told to
    
    Parameters
    ----------
    neighborList : 
        the neighbor list object
    pot :
        the potential object
    """

    def __init__(self, neighborList, pot):
        self.neighborList = neighborList
        self.pot = pot

    def buildList(self, coords):
        """
        instruct the neighbor list object to rebuild it's list
        """
        self.list = self.neighborList.buildList(coords)

    def getEnergy(self, coords):
        return self.pot.getEnergyList(coords, self.list)

    def getEnergyGradient(self, coords):
        return self.pot.getEnergyGradientList(coords, self.list)


class NeighborListPotentialMulti(basepot):
    """
    A wrapper for multiple NeighborListPotentialBuild
    
    This will wrap multiple instances of NeighborListPotentialBuild.  This class
    will the check the coords and control when the neighbor lists are rebuilt.
    
    Parameters:
    -----------
    potentials :
        a list of neighbor list potential objects.  Each potential must have attribute
        
            potential.buildList(coords)
    natoms : 
        number of atoms
    rcut : 
        the cutoff distance for the potential
    rskin : 
        the skin distance.  atoms are listed as 
        neighbors if they are closer then rlist = rcut + rskin.  A larger rskin
        will have more interaction pairs, but will need to be updated less
        frequently
    boxl : 
        if not None, then the system is in a periodic box of size boxl
        
    """

    def __init__(self, potentials, natoms, rcut, rskin=0.5, boxl=None):
        self.potentials = potentials
        self.oldcoords = np.zeros([natoms, 3])
        self.rcut = rcut
        self.rskin = rskin
        self.redo_displacement = self.rskin / 2.
        self.rlist = self.rcut + self.rskin
        self.rlist2 = self.rlist ** 2

        if boxl is None:
            self.periodic = False
            self.boxl = 1.
        else:
            self.periodic = True
            self.boxl = boxl

        self.buildcount = 0
        self.count = 0

    def needNewList(self, coords):
        coords = np.reshape(coords, [-1, 3])
        if self.periodic:
            #only check periodic boundary conditions for the atoms that fail the normal test
            indices = np.where(((coords - self.oldcoords) ** 2).sum(1) > self.redo_displacement ** 2)[0]
            if len(indices) == 0:
                return False
            dr = coords[indices, :] - self.oldcoords[indices, :]
            dr -= self.boxl * np.round(dr / self.boxl)
            return np.any((dr ** 2).sum(1) > self.redo_displacement ** 2)
        else:
            return np.any(((coords - self.oldcoords) ** 2).sum(1) > self.redo_displacement ** 2)

    def update(self, coords):
        self.count += 1
        if self.needNewList(coords):
            self.buildcount += 1
            self.oldcoords = np.copy(coords).reshape([-1, 3])
            for pot in self.potentials:
                pot.buildList(coords)

    def getEnergy(self, coords):
        self.update(coords)
        E = 0.
        for pot in self.potentials:
            E += pot.getEnergy(coords)
        return E

    def getEnergyGradient(self, coords):
        self.update(coords)
        Etot = 0.
        gradtot = np.zeros(np.shape(coords))
        for pot in self.potentials:
            E, grad = pot.getEnergyGradient(coords)
            Etot += E
            gradtot += grad
        return Etot, gradtot


class MultiComponentSystem(basepot):
    """
    a potential wrapper for multiple potentials
    
    Parameters
    ----------
    potentials :
        a list of potential objects
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

#def makeBLJNeighborListPot(natoms, ntypeA = None, rcut=2.5, boxl=None):
#    """
#    recreate the binary lj with atom typea A,B from 3 interaction lists AA, BB, AB
#    """
#    print "making BLJ neighborlist potential", natoms, ntypeA, rcut, boxl
#    #rcut = 2.5
#    #natoms = 40
#    if ntypeA is None:
#        ntypeA = int(natoms * 0.8)
#    Alist = range(ntypeA)
#    Blist = range(ntypeA, natoms)
#    
#    
#    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut)
#
#        
#    ljAA = LJ(eps=blj.AA.eps, sig=blj.AA.sig, rcut=rcut*blj.AA.sig, boxl=boxl)
#    ljBB = LJ(eps=blj.BB.eps, sig=blj.BB.sig, rcut=rcut*blj.BB.sig, boxl=boxl)
#    ljAB = LJ(eps=blj.AB.eps, sig=blj.AB.sig, rcut=rcut*blj.AB.sig, boxl=boxl)
#    
#    nlAA = NeighborListSubset(natoms, rcut, Alist, boxl=boxl )
#    nlBB = NeighborListSubset(natoms, rcut, Blist, boxl=boxl )
#    nlAB = NeighborListSubset(natoms, rcut, Alist, Blist, boxl=boxl)
#    
#    potlist = [ 
#               NeighborListPotential(nlAA, ljAA),
#               NeighborListPotential(nlBB, ljBB),
#               NeighborListPotential(nlAB, ljAB)
#                ]
#    mcpot = MultiComponentSystem(potlist)
#    return mcpot



#def test(natoms = 40, boxl=None):
#    import pele.potentials.ljpshiftfast as ljpshift
#    from pele.optimize import mylbfgs
#    ntypeA = int(natoms*0.8)
#    rcut = 2.5
#    coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2
#    
#    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)
#
#    pot = makeBLJNeighborListPot(natoms, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
#    
#
#    
#    eblj = blj.getEnergy(coords)
#    print "blj energy", eblj
#    
#    epot = pot.getEnergy(coords)
#    print "mcpot energy", epot
#    
#    print "difference", (epot - eblj)/eblj
#    
#    ret1 = mylbfgs(coords, blj, iprint=-11)
#    np.savetxt("out.coords", ret1.coords)
#    print "energy from quench1", ret1.energy
#    ret2 = mylbfgs(coords, pot, iprint=-1)
#    print "energy from quench2", ret2.energy
#    
#    print "ret1 evaluated in both potentials", pot.getEnergy(ret1.coords), blj.getEnergy(ret1.coords)
#    print "ret2 evaluated in both potentials", pot.getEnergy(ret2.coords), blj.getEnergy(ret2.coords)
#    
#    coords = ret1.coords
#    e1, g1 = blj.getEnergyGradient(coords)
#    e2, g2 = pot.getEnergyGradient(coords)
#    print "energy difference from getEnergyGradient", (e2 - e1)
#    print "largest gradient difference", np.max(np.abs(g2-g1))
#    print "rms gradients", np.linalg.norm(g1)/np.sqrt(len(g1)), np.linalg.norm(g2)/np.sqrt(len(g1))
#
#    
#    for subpot in pot.potentials:
#        nl = subpot.neighborList
#        print "number of times neighbor list was remade", nl.buildcount, "out of", nl.count
#    
#    if False:
#        try: 
#            import pele.utils.pymolwrapper as pym
#            pym.start()
#            pym.draw_spheres(np.reshape(coords,[-1,3]), "A", 1)
#            pym.draw_spheres(np.reshape(ret1.coords,[-1,3]), "A", 2)
#            pym.draw_spheres(np.reshape(ret2.coords,[-1,3]), "A", 3)
#        except ImportError:
#            print "Could not draw using pymol, skipping this step" 
#    
#def test2():
#    import pele.potentials.ljpshiftfast as ljpshiftfast
#    import pele.potentials.ljpshift as ljpshift
#    from pele.optimize import mylbfgs
#    fname = "/scratch/scratch2/js850/library/cluster/spherical/1620/PTMC/q4/oneatom/cavity200-8/ts/coords1.quench"
#    fname = "/scratch/scratch2/js850/library/cluster/spherical/1620/PTMC/q4/oneatom/cavity200-8/ts/test.coords"
#    #fname = "out.coords"
#    if False:
#        coords = np.array(np.loadtxt(fname))
#        coords = coords.reshape(-1)
#        boxl = 11.05209
#    else:
#        natoms = 200
#        coords = np.random.uniform(-1,1,natoms*3)*(natoms)**(1./3)/2
#        print "max, min coords", coords.max(), coords.min()
#        boxl = 5
#
#    natoms = len(coords) /3
#    ntypeA = int(natoms*0.8)
#    rcut = 2.5
#    print "natoms", natoms, "ntypea", ntypeA
#    
#    blj = ljpshift.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)
#    bljfast = ljpshiftfast.LJpshift(natoms, ntypeA, rcut=rcut, boxl=boxl)
#
#    pot = makeBLJNeighborListPot(natoms, ntypeA=ntypeA, rcut=rcut, boxl=boxl)
#    
#    eblj = blj.getEnergy(coords)
#    print "blj energy", eblj
#    
#    epot = pot.getEnergy(coords)
#    print "mcpot energy", epot
#    
#    print "energy difference", (epot - eblj)
#    
#    e1, g1 = blj.getEnergyGradient(coords)
#    e2, g2 = pot.getEnergyGradient(coords)
#    print "energy difference from getEnergyGradient", (e2 - e1)
#    print "largest gradient difference", np.max(np.abs(g2-g1))
#    print "rms gradients", np.linalg.norm(g1)/np.sqrt(len(g1)), np.linalg.norm(g2)/np.sqrt(len(g1))
#    
#    if False:
#        print "quenching"
#        ret1 = mylbfgs(coords, blj, iprint=-11)
#        np.savetxt("out.coords", ret1.coords)
#        print "energy from quench1", ret1.energy
#        ret2 = mylbfgs(coords, pot, iprint=-1)
#        print "energy from quench2", ret2.energy
#        print "max, min quenched coords", coords.max(), coords.min()
#
#
#        print "ret1 evaluated in both potentials", pot.getEnergy(ret1.coords), blj.getEnergy(ret1.coords)
#        print "ret2 evaluated in both potentials", pot.getEnergy(ret2.coords), blj.getEnergy(ret2.coords)
#    elif True:
#        print "quenching"
#        ret2 = mylbfgs(coords, pot, iprint=-1)
#        print "energy from quench2", ret2.energy
#        print "max, min quenched coords", ret2.coords.max(), ret2.coords.min()
#
#        print "ret2 evaluated in both potentials", pot.getEnergy(ret2.coords), blj.getEnergy(ret2.coords)
#        print "and in blj fast                  ", bljfast.getEnergy(ret2.coords)
#
#        
#
#
#    
#if __name__ == "__main__":
#    #test2()
#    test(natoms=100, boxl=None)
    
       
        

