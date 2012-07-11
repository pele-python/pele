import lattice
import vec3
from pygmin import rotations
import numpy as np

def compareTransformed(coords1, coords2, x1, x2, M):
    # the lattice matrix
    ml1 = lattice.lowerTriangular(coords1.lattice)
    ml2 = lattice.lowerTriangular(coords2.lattice)
    
    # the inerse lattice matrix
    iml2 = vec3.invert3x3(ml2)

    for i in xrange(coords1.nrigid):
        ptest = coords1.rotRigid[i]
        ptrans = rotations.mx2aa(rotations.np.dot(M, rotations.aa2mx(ptest)))
        for d in [np.zeros(3),[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]:
            #print "Testing displace",d
            xtest = coords1.posRigid[i] + np.array(d) - x1
            xtrans = np.dot(iml2, np.dot(M, np.dot(ml1, xtest)))
            #print xtrans
            match = False
            for j in xrange(coords1.nrigid):
                dx = coords2.posRigid[j] - xtrans - x2
                dx = dx - np.floor(dx)
                dx[dx>0.8]-=1.0
                dp = rotations.mx2aa(np.dot(vec3.invert3x3(rotations.aa2mx(ptrans)),rotations.aa2mx(coords2.rotRigid[j])))
                #print i,j,np.linalg.norm(dx), np.linalg.norm(dp)
                # print dx, dp
                match = np.linalg.norm(dx) < 1e-2 \
                    and np.linalg.norm(dp) < 1e-3
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
            match = compareTransformed(coords1, coords2, x1, x2, M)
            if(match):
                #print "Structures match"
                return True
    print "no match"
    #import dmagmin_ as GMIN
    #GMIN.writeCIF("1.cif", coords1.coords)
    #GMIN.writeCIF("2.cif", coords2.coords)
    #import pickle
    #pickle.dump(coords1, open("1.dat", "w"))
    #pickle.dump(coords2, open("2.dat", "w"))
    #exit()
    return False
            