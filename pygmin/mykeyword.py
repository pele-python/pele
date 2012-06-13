import numpy as np
import basinhopping as bh

class myKeywordClass():
    def __init__(self):
        self.setDefaults()

    def setDefaults(self):
        self.pot_type = "lj"
        self.sig = 1. 
        self.eps = 1.
        self.periodic = False
        self.boxl = None
        self.stepsize = 0.3
        self.nmcsteps = 100
        self.temperature = 1.0
        self.accrat = 0.5
        self.accrat_frq = 50
        self.nsave = 1
        self.ediff = 0.002

    def getBasinHopping(self, coords, natoms):
        #########################################################################
        #initialize potential, etc
        #########################################################################
        if self.pot_type == "lj":
            import potentials.lj as lj
            potential = lj.LJ(self.eps, self.sig, self.boxl)
        if self.pot_type == "binary":
            usefortran = True
            if usefortran:
                import potentials.ljpshiftfast as ljpshift
            else:
                import potentials.ljpshift as ljpshift
            potential = ljpshift.LJpshift( natoms, self.ntypeA, self.boxl, self.cutoff, self.epsBB, self.sigBB, self.epsAB, self.sigAB)

        event_after_step = []

        #initialize minima saving routine
        import storage.savenlowest as saveit
        self.savelowest = saveit.SaveN( nsave = self.nsave, accuracy = self.ediff ) #class to save the lowest energy structures


        #initialize step taking routine using adaptive step size class
        import take_step.random_displacement as random_displacement
        self.takeStep = random_displacement.takeStep( ) #class to implement the take step routine
        #initialize adaptive step size routine
        self.takeStep.useAdaptiveStep(stepsize=self.stepsize, acc_ratio = self.accrat, freq = self.accrat_frq )

        #classes to impiment acceptence criterion
        acceptTests = []
        #add further tests here, e.g    acceptTests.append( cold_fusion_test )

        #add optional events, e.g. dump coords
        #event_after_step.append(  myDumpRoutine )

        #initialize basing hopping class and return it
        opt = bh.BasinHopping(coords, potential, \
                takeStep = self.takeStep, \
                storage = self.savelowest.insert,  \
                event_after_step=event_after_step, \
                acceptTests=acceptTests, \
                temperature=self.temperature \
                )
        return opt

    def readKeywords(self, fin):
        for line in fin:
            self.addKeyword(line)

    def addKeyword(self, line):
            line = line.lower() #put in lower case
            words = line.split()
            if words[0] == "comment":
                return
            elif line[0] == "#":
                return

            elif words[0] == "lj":
                """lennard jones potential"""
                self.pot_type = "lj"
                self.eps = float(words[1])
                self.sig = float(words[2])

            elif words[0] == "blj":
                """binary lennard jones potential"""
                self.pot_type = "blj"
                self.ntypeA = int(words[1])
                self.epsAB = float(words[2])
                self.epsBB = float(words[3])
                self.sigAB = float(words[4])
                self.sigBB = float(words[5])

            elif words[0] == "binary":
                """binary lennard jones potential with cutoff"""
                self.pot_type = "binary"
                self.ntypeA = int(words[1])
                self.epsAB = float(words[2])
                self.epsBB = float(words[3])
                self.sigAB = float(words[4])
                self.sigBB = float(words[5])

            elif words[0] == "shiftcut":
                """cutoff distance for the potential"""
                self.cutoff = float(words[1])

            elif words[0] == "periodic":
                """box size for periodic systems (cubic only at the moment)"""
                self.periodic = True
                self.boxl = float(words[1])

            elif words[0] == "step":
                """initial step size for the monte carlo step"""
                self.stepsize = float(words[1])

            elif words[0] == "steps":
                """number of monte carlo steps"""
                self.nmcsteps = int(words[1])

            elif words[0] == "accrat" or words[0] == "acceptratio":
                """target acceptance ratio"""
                self.accrat = float(words[1])

            elif words[0] == "changeaccept" or words[0] == "acceptratio":
                """the interval at which the acceptance ratio is modified"""
                self.accrat_frq = float(words[1])

            elif words[0] == "temperature":
                """temperature for the monte carlo acceptance criterion"""
                self.temperature = float(words[1])

            elif words[0] == "save":
                """number of lowest minima to save"""
                self.nsave = int(words[1])

            elif words[0] == "ediff":
                """energy difference used to distinguish different minima in the save lowest routine"""
                self.ediff = float(words[1])


            else:
                print "keyword ", words[0], " not implimented"



