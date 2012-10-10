__all__ = ["printAtomsXYZ", "PrintEvent"]

def printAtomsXYZ(fout, coords, line2="", atom_type=["LA"]):
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
#    if isinstance(atom_type, str):
#        atomtype = [atom_type]
    natomtypes = len(atom_type)
    fout.write( str(natoms) + "\n")
    fout.write( str(line2) + "\n")
    for i in xrange(natoms):
        fout.write( atom_type[i % natomtypes] +" "+ str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 



class PrintEvent(object):
    def __init__(self, fout, frq = 1, atom_type=["LA"]):
        """
        for passing a print event to minimizers and other objects
        """
        if isinstance(fout, str): 
            fname = fout      
            self.fout = open(fname, "w")
        else:
            self.fout = fout
        self.atom_type = atom_type
        self.frq = frq
        self.count = 0
    
    def event(self, e, coords, accepted):
        self.count += 1
        if self.count % self.frq == 0:
            printAtomsXYZ(self.fout, coords, str(self.count) + " " + str(e))
    
    def __call__(self, e, coords, accepted):
        self.event(e, coords, accepted)
        