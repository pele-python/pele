"""
.. currentmodule:: pygmin.utils.crystals

.. autosummary::
    :toctree: generated/
    
    compareTransformed
    compareStructures
    has_overlap
"""

import lattice
import vec3
from pygmin.utils import rotations
import numpy as np
import vec3

__all__ = ["compareTransformed", "compareStructures", "has_overlap"]

tol_rot = 3.1415 / 180. # standard tolerance is 1 deg
tol_shift = 0.01 # standard tolerance is 0.1 in absolute coordinates

GMIN=None

def compareTransformed(coords1, coords2, x1, x2, M):
    # the lattice matrix
    ml1 = lattice.lowerTriangular(coords1.lattice)
    ml2 = lattice.lowerTriangular(coords2.lattice)
    
    # the inerse lattice matrix
    iml2 = vec3.invert3x3(ml2)
    for i in xrange(coords1.nrigid):
        ptest = coords1.rotRigid[i]
        ptrans = rotations.mx2aa(rotations.np.dot(M, rotations.aa2mx(ptest)))
        for d in [np.zeros(3),[0.,1.,0.],[0.,1.,0.],[0.,0.,1.]]:
            xtest = coords1.posRigid[i] + np.array(d) - x1
            xtrans = np.dot(iml2, np.dot(M, np.dot(ml1, xtest)))
     
            match = False
            for j in xrange(coords1.nrigid):
                dx = coords2.posRigid[j] - xtrans - x2
                dx = dx - np.floor(dx)
                dx[dx>0.8]-=1.0
                dp = rotations.mx2aa(np.dot(vec3.invert3x3(rotations.aa2mx(ptrans)),rotations.aa2mx(coords2.rotRigid[j])))
                match = np.linalg.norm(np.dot(ml2, dx)) < tol_shift \
                    and np.linalg.norm(dp) < tol_rot
                if(match):
                    break
            if(not match):
                return False
    return True

def compareStructures(coords1, coords2):    
    for refmol1 in xrange(coords1.nrigid):
        for refmol2 in xrange(refmol1, coords1.nrigid):
            x1 = coords1.posRigid[refmol1]
            x2 = coords2.posRigid[refmol2]
            pref1 = coords1.rotRigid[refmol1]
            pref2 = coords2.rotRigid[refmol2]
            #print "pref",pref1,pref2
            M = vec3.invert3x3(rotations.aa2mx(pref1))
            M = np.dot(rotations.aa2mx(pref2), M)
            # print x1,x2
            ptest = pref1
            ptrans = rotations.mx2aa(rotations.np.dot(M, rotations.aa2mx(ptest)))
            #print "---------------- start try",refmol1,refmol2
            match1 = compareTransformed(coords1, coords2, x1, x2, M)
            match2 = compareTransformed(coords2, coords1, x2, x1, vec3.invert3x3(M))
            if(match1 and match2):
                #print "Structures match"
                return True
    #import dmagmin_ as GMIN
    #GMIN.writeCIF("1.cif", coords1.coords)
    #GMIN.writeCIF("2.cif", coords2.coords)
    #import pickle
    #pickle.dump(coords1, open("1.dat", "w"))
    #pickle.dump(coords2, open("2.dat", "w"))
    #exit()
    return False

def has_overlap(coords, cutoff=1.0):
    ml = lattice.lowerTriangular(coords[-6:])
    natoms = coords.size/3-2
    atoms = coords.reshape(coords.size/3,3)[:-2]
    # mi = vec3.invert3x3(ml)
    
    mindist = 1000.

    a = ml[:,0]
    b = ml[:,1]
    c = ml[:,2]

    for i in xrange(natoms):
        for j in xrange(i+1,natoms):
            #d = np.dot(mi, atoms[i] - atoms[j])
            #d -= np.floor(d)
            #d[d>0.5] -= 1.
            #r_ij = np.linalg.norm(np.dot(ml, d))
            
            r_i = atoms[i]
            r_j = atoms[j]
            r_tp = r_j - r_i;
            r_dp = r_tp - c*int(r_tp[2]/c[2])
            r_sp = r_dp - b*int(r_dp[1]/b[1])
            r_ij = r_sp - a*int(r_sp[0]/a[0])

            r = np.linalg.norm(r_ij)
            if(r < cutoff):
                # print "overlapping distance is", r
                return True
            mindist = min(mindist, r)
                
    #print "the minimum distance ist", mindist
    return False
    
    