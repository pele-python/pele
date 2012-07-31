import numpy as np
from mindistutils import CoMToOrigin, alignRotation

def minDist(X1, X2):
    """
    Minimize the distance between two clusters.  The following symmetries will be accounted for
    
    Translational symmetry

    Global rotational symmetry
    """
    #alignCoM(X1, X2)
    X1 = CoMToOrigin(X1)
    X2 = CoMToOrigin(X2)

    #align rotation degrees of freedom
    dist, X2 = alignRotation(X1, X2)
    return dist, X1, X2


def main():
    natoms = 5
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

    dist, X1, X2 = minDist(X1,X2)
    distfinal = np.linalg.norm(X1-X2)
    print "dist from eigenvalue", dist
    print "distfinal", distfinal

    import printing.print_atoms_xyz as printxyz
    with open("out.xyz", "w") as fout:
        CoMToOrigin(X1i)
        CoMToOrigin(X2i)
        printxyz.printAtomsXYZ(fout, X1i )
        printxyz.printAtomsXYZ(fout, X2i )
        printxyz.printAtomsXYZ(fout, X1 )
        printxyz.printAtomsXYZ(fout, X2 )

if __name__ == "__main__":
    main()
