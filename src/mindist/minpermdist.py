import numpy as np
import itertools
import rotations as rot
import mindist
import math


def permuteArray(Xold, perm):
    #don't modify Xold
    Xnew = np.copy(Xold)
    for (iold, inew) in enumerate(perm):
        #print iold, "->", inew
        Xnew[inew*3:inew*3+3] = Xold[iold*3:iold*3+3]

    return Xnew

def findBestPermutation( X1, X2):
    """
    For a given set of positions X1 and X2, find the best permutation of the
    atoms in X2.

    Use an implimentation of the Hungarian Algorithm in the Python package
    index (PyPi) called munkres (another name for the algorithm).  The
    hungarian algorithm time scales as O(n^3), much faster than the O(n!) from
    looping through all permutations.

    http://en.wikipedia.org/wiki/Hungarian_algorithm
    http://pypi.python.org/pypi/munkres/1.0.5.2
    """
    try:
        from munkres import Munkres
    except:
        print "munkres package not installed, skipping Hungarian algorithm"
        dist = np.linalg.norm( X1 - X2 )
        return dist, X1, X2

    nsites = len(X1) / 3


    #########################################
    # create the cost matrix
    #########################################
    cost = np.zeros( [nsites,nsites], np.float64)
    for i in range(nsites):
        for j in range(nsites):
            R2 = np.sum( (X1[i*3:i*3+3] - X2[j*3:j*3+3])**2 )
            #R2old = np.linalg.norm( X1[i*3:i*3+3] - X2[j*3:j*3+3] )**2
            #if abs(R2 - R2old) > 1e-6:
                #print "i'm going crazy", R2, R2old
            cost[j,i] = R2

    #convert cost matrix to a form used by munkres
    matrix = cost.tolist()

    #########################################
    # run the hungarian algorithm
    #########################################
    m = Munkres()
    newind = m.compute(matrix)

    #########################################
    # apply the permutation
    #########################################
    costnew = 0.;
    X2old = np.copy(X2)
    for (iold, inew) in newind:
        costnew    += cost[iold,inew]
        if iold != inew:
            #print iold, "->", inew, "matrix %10.4f, %10.4f" % (matrix[iold][inew], matrix[inew][iold])
            #for i in [iold, inew]:
                #for j in [iold, inew]:
                    #r = np.linalg.norm( X1[i*3:i*3+3] - X2old[j*3:j*3+3] )
                    #print "    %4d %4d %10.4f, %10.4f" % (i, j, r, r**2 ), cost[j, i], matrix[i][j]
            X2[inew*3:inew*3+3] = X2old[iold*3:iold*3+3]

    #costold = sum( [matrix[i][i] for i in range(nsites)] )
    #print "costold    ", costold, np.sqrt(costold)
    #print "costnew    ", costnew, np.sqrt(costnew)

    dist = np.sqrt(costnew)
    return dist, X1, X2

def minPermDistLong(X1, X2, max_permutations = 10000):
    """
    Minimize the distance between two clusters.  The following symmetries will be accounted for
    
    Translational symmetry

    Global rotational symmetry

    Permutational symmetry
    """
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

        dist = mindist.minDist(X1, X2)

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

    import printing.print_atoms_xyz as printxyz
    with open("out.xyz", "w") as fout:
        mindist.CoMToOrigin(X1i)
        mindist.CoMToOrigin(X2i)
        printxyz.printAtomsXYZ(fout, X1i )
        printxyz.printAtomsXYZ(fout, X2i )
        printxyz.printAtomsXYZ(fout, X1 )
        printxyz.printAtomsXYZ(fout, X2 )

if __name__ == "__main__":
    main()
