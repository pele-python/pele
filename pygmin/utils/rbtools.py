class CoordsAdapter():
    def __init__(self, nrigid=0, natoms=0, nlattice=0, coords=None):
        self.nrigid = nrigid
        self.natoms = natoms
        self.nlattice = nlattice
        if(coords!=None):
            self.updateCoords(coords)
        
    def updateCoords(self, coords):
        natoms = self.natoms
        nrigid = self.nrigid
        self.coords = coords
        
        self.posAtoms=None
        self.posRigid=None
        self.rotRigid=None
        self.lattce=None
        
        if natoms > 0:
            self.posAtoms = self.coords[6*nrigid:6*nrigid+3*natoms].reshape(natoms, 3)
        
        if nrigid > 0:
            self.posRigid = self.coords[0:3*nrigid].reshape(nrigid, 3)
            self.rotRigid = self.coords[3*nrigid:6*nrigid].reshape(nrigid, 3)
    
        if self.nlattice > 0:
            self.lattice = self.coords[-self.nlattice:]
