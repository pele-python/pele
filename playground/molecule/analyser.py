from playground.molecule.molecule import Molecule, Protein

class Analyser(object):
    ''' Class provides a suite of analysis routines for molecular objects '''
    def __init__(self, *args, **kwargs):
        ''' Initialise an analyser object '''
        print "Molecular Analyser Initialised"
        
    def __analyse__(self, molecule):
        ''' This function analyses a molecule object'''
        print "Analysing Molecule"
        print molecule.name


class ProteinAnalyser(Analyser):
    ''' Class provides a suite of analysis routines for protein molecules '''
    def __init__(self, *args, **kwargs):
        ''' Initialise an ProteinAnalyser object '''
        print "Protein Analyser Initialised"
        
    def __analyse__(self, protein):
        ''' This function analyses a molecule object'''
        print "Analysing Protein"
        print protein.name
    
