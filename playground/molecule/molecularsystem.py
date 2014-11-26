class MolecularSystem(object):
    ''' A container for molecules which are defined as single component graphs'''
    
    def __init__(self, molecules):
        ''' Initialise the molecular system with a list of molecule objects '''
        self.molecules = molecules
        self.number_molecules = len(molecules)
    