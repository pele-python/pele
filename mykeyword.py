
class myKeywordClass():
    def __init__(self):
        a=1

    def setDefaults(self):
        self.potential = "lj"
        self.sig = 1. 
        self.eps = 1.
        self.periodic = False
        self.boxl = None
        self.stepsize = 0.3
        self.nmcsteps = 100
        self.temperature = 1.0
        self.accrat = 0.5
        self.accrat_frq = 50

    def readKeywords(self, fin):
        for line in fin:
            line = line.lower() #put in lower case
            words = line.split()
            if words[0] == "comment":
                continue
            elif line[0] == "#":
                continue

            elif words[0] == "lj":
                self.potential = "lj"
                self.eps = float(words[1])
                self.sig = float(words[2])

            elif words[0] == "ljcpp":
                self.potential = "ljcpp"

            elif words[0] == "blj":
                self.potential = "blj"
                self.ntypeA = int(words[1])
                self.epsAB = float(words[2])
                self.epsBB = float(words[3])
                self.sigAB = float(words[4])
                self.sigBB = float(words[5])

            elif words[0] == "binary":
                self.potential = "binary"
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



