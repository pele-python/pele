def printAtomsXYZ(fout, coords, line2="", atom_type="LA"):
    """
    print atomic coordinate in vmd xyz format:

    natoms
    line2 ...could be anything, e.g. box lengths...
    atom_type x[0] y[0] z[0]
    atom_type x[1] y[1] z[1]
    atom_type x[2] y[2] z[2]
    ...
    """
    natoms = len(coords)/3
    fout.write( str(natoms) + "\n")
    fout.write( str(line2) + "\n")
    for i in xrange(natoms):
        fout.write( atom_type +" "+ str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 


