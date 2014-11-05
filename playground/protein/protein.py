import molecule


class Protein(Molecule):
    ''' A representation of a protein that extends the Molecule base class.
        The basic graph that describe the molecule in the base class is analysed to determine:
        
        Chains
        Residues
        Side groups
        Backbone '''

    def init(self,l_input_filename):
        # Call the initialisation routine of the base class
        super(self.__class__, self).__init__(l_input_filename)
        
        # Analyse graph to determine chains
        self.chains=self.identify_chains()

        #self.residues=self.
        
    
    def identifiy_chains(self):
        
        return chains
    
