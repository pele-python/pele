import numpy as np
import itertools
from pygmin.mindist import minDist
import math
from pygmin.mindist import findBestPermutation, CoMToOrigin
from pygmin.mindist.mindistutils import permuteArray

__all__ = ["minPermDistLong"]

def minPermDistLong(X1, X2, max_permutations = 10000):
    """
    Minimize the distance between two clusters.  The following symmetries will be accounted for
    
    Translational symmetry

    Global rotational symmetry

    Permutational symmetry
    
    This routine is deterministic, but ludicrously slow.  Use the minPermDistStochastic instead
    """
    print "This routine is deterministic, but ludicrously slow.  Use the minPermDistStochastic instead"

    X2in = np.copy(X2)

    X2min = np.copy(X2)
    dmin = np.linalg.norm( X2-X1 )

    nsites = len(X1) / 3
    #this is a really dumb way to do this
    nperm = math.factorial(nsites)
    print nperm, "permutations in total. Stopping at ", max_permutations
    count = 0
    for perm in itertools.permutations(range(nsites)):
        #print perm
        X2 = permuteArray( X2in, perm)

        dist, X1, X2 = minDist(X1, X2)

        if dist < dmin:
            dmin = dist
            X2min = np.copy(X2)
            #print dist, np.linalg.norm(X1-X2)
            if dmin < 1e-6:
                print "found exact match (to accuracy 1e-6)"
                break


        count += 1
        if count % 1000 == 0:
            print "finished", count, "permutations out of", min(nperm, max_permutations), " mindist is", dmin, "last try was", dist



        if count == max_permutations:
            print "reached maximum number of perterbations (", max_permutations, ")."  

            use_hungarian = True
            if use_hungarian:
                print "will now fix alignment and calculate the best permutation using the Hungarian algorithm"
                #print "dist before hungarian algorithm", dmin
                ret = findBestPermutation( X1, X2min )
                if ret != None:
                    dmin, X1, X2min = findBestPermutation( X1, X2min )
                #print "dist after hungarian algorithm", dmin


            break


    return dmin, X1, X2min


def main():
    natoms = 9
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

    (dist, X1, X2) = minPermDistLong(X1,X2)
    distfinal = np.linalg.norm(X1-X2)
    print "dist returned    ", dist
    print "dist from coords ", distfinal

    import pygmin.printing.print_atoms_xyz as printxyz
    with open("out.xyz", "w") as fout:
        CoMToOrigin(X1i)
        CoMToOrigin(X2i)
        printxyz.printAtomsXYZ(fout, X1i )
        printxyz.printAtomsXYZ(fout, X2i )
        printxyz.printAtomsXYZ(fout, X1 )
        printxyz.printAtomsXYZ(fout, X2 )

if __name__ == "__main__":
    main()
