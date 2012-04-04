import numpy as np
import basinhopping as bh
import saveit
import adaptive_step
import take_step

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

    def getBasinHopping(self, coords, natoms):
        #########################################################################
        #initialize potential, etc
        #########################################################################
        if self.pot_type == "lj":
            import potentials.lj as lj
            potential = lj.LJ(self.eps, self.sig, self.boxl)
        if self.pot_type == "ljcpp":
            import potentials.ljcpp as ljcpp
            potential = ljcpp.LJ()
        if self.pot_type == "binary":
            usefortran = True
            if usefortran:
                import potentials.ljpshiftfast as ljpshift
            else:
                import potentials.ljpshift as ljpshift
            potential = ljpshift.LJpshift( natoms, self.ntypeA, self.boxl, self.cutoff, self.epsBB, self.sigBB, self.epsAB, self.sigAB)

        #initialize minima saving routine
        self.savelowest = saveit.saveit() #class to save the lowest energy structure

        #initialize adaptive step size routine
        self.manstep = adaptive_step.manageStepSize (self.stepsize, self.accrat, self.accrat_frq) #class to do step size adjustment
        event_after_step = [self.manstep.insertStepWrapper]


        #initialize step taking routine using adaptive step size class
        takeStep = take_step.takeStep( RNG = np.random.rand, getStep = self.manstep.getStepSize ) #class to impliment the take step routine

        #classes to impiment acceptence criterion
        acceptTests = []
        #add further tests here, e.g    acceptTests.append( cold_fusion_test )

        #add optional events, e.g. dump coords
        #event_after_step.append(  myDumpRoutine )

        #initialize basing hopping class and return it
        opt = bh.BasinHopping(coords, potential, takeStep, storage = self.savelowest.insert,  \
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
                self.pot_type = "lj"
                self.eps = float(words[1])
                self.sig = float(words[2])

            elif words[0] == "ljcpp":
                self.pot_type = "ljcpp"

            elif words[0] == "blj":
                self.pot_type = "blj"
                self.ntypeA = int(words[1])
                self.epsAB = float(words[2])
                self.epsBB = float(words[3])
                self.sigAB = float(words[4])
                self.sigBB = float(words[5])

            elif words[0] == "binary":
                self.pot_type = "binary"
                self.ntypeA = int(words[1])
                self.epsAB = float(words[2])
                self.epsBB = float(words[3])
                self.sigAB = float(words[4])
                self.sigBB = float(words[5])

            elif words[0] == "shiftcut":
                self.cutoff = float(words[1])

            elif words[0] == "periodic":
                self.periodic = True
                self.boxl = float(words[1])

            elif words[0] == "step":
                self.stepsize = float(words[1])

            elif words[0] == "steps":
                self.nmcsteps = int(words[1])

            elif words[0] == "accrat" or words[0] == "acceptratio":
                self.accrat = float(words[1])

            elif words[0] == "changeaccept" or words[0] == "acceptratio":
                self.accrat_frq = float(words[1])

            elif words[0] == "temperature":
                self.temperature = float(words[1])


            else:
                print "keyword ", words[0], " not implimented"



