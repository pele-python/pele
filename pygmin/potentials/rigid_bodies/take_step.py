import numpy as np
import pygmin.utils.rotations as rot


class RBTakeStep:
    """
    a take step routine for rigid bodies
    
    warning: adaptive step size not implemented
    """
    def __init__(self, RNG = np.random.rand, \
                 getTStep=None, Tstepsize=0.3, \
                 getOStep=None, Ostepsize=np.pi/4. \
                 ):
        self.RNG = RNG #random number generator
        if getTStep == None: 
            self.getTStep = lambda : Tstepsize
        else:
            print "warning: adaptive step size not implemented"
            self.getTStep = getTStep
    
        if getOStep == None: 
            self.getOStep = lambda : Ostepsize
        else:
            print "warning: adaptive step size not implemented"
            self.getOStep = getOStep
            
        self.stepnum = 0
        #the user can modify steporder
        self.steporder = [self.translation_step, self.orientational_step]
        self.nstep_types = len(self.steporder)
    
    def takeStep(self, coords, **kwargs):
        """
        take steps in the order defined by self.steporder
        """
        this_step = self.steporder[ self.stepnum % self.nstep_types ]
        this_step(coords)
        self.stepnum += 1
    
    def translation_step(self, coords):
        """
        take a random translational step.

        note: this step is rotationally invariant, which is different than the normal take_step routine
        note: this method also favors shorter steps than the normal take_step routine
        """
        nmol = len(coords) / 2 / 3
        stepsize = self.getTStep()
        for j in range(nmol):
            #get a random unit vector
            #this is a dumb way to do it
            v = rot.random_aa()
            #make the step length random and uniform in [0,stepsize]
            v *= self.RNG() * stepsize / np.linalg.norm(v)
            #print "rand ", rand
            coords[3*j : 3*j + 3] += v
            
    def orientational_step(self, coords):
        """
        take a random orientational step
        """
        nmol = len(coords) / 2 / 3
        maxtheta = self.getOStep()
        for j in range(nmol):
            rot.takestep_aa( coords[3*nmol + 3*j : 3*nmol + 3*j + 3], maxtheta )

    def updateStep(self, accepted, **kwargs):
        pass


def test_takestep(nmol = 4):
    #get an initial set of coordinates
    import copy
    
    nsites = nmol*3
    comcoords = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    aacoords = np.array( [copy.copy(rot.random_aa()) for i in range(nmol)] )
    aacoords = aacoords.reshape(3*nmol)
    coords = np.zeros(2*3*nmol, np.float64)
    coords[0:3*nmol] = comcoords[:]
    coords[3*nmol:2*3*nmol] = aacoords[:]
    print "lencoords, len aacoords", len (coords), len(aacoords), len(comcoords)

    takestep = RBTakeStep()
    print coords
    oldcoords = coords.copy()
    takestep.translation_step(coords)
    print coords
    print coords - oldcoords
    
    print ""
    print "taking orientational step"
    print ""
    oldcoords = coords.copy()
    takestep.orientational_step(coords)
    print coords
    print coords - oldcoords
    
    print ""
    print "taking one of each step"
    print ""
    oldcoords = coords.copy()
    takestep.takeStep(coords)
    takestep.takeStep(coords)
    print coords
    print coords - oldcoords



if __name__ == "__main__":
    test_takestep()

