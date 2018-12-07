from __future__ import print_function
from pele.takestep import TakestepInterface

__all__ = ["AdaptiveStepsizeTemperature"]


class AdaptiveStepsizeTemperature(TakestepInterface):
    """
    adjust both the stepsize and the temperature adaptively
    
    Parameters
    ----------
    stepclass : 
        the step taking class
    target_new_min_prob : 
        the target probability for a step ending up in a new minimum.
        Used to adjust the stepsize.
    target_new_min_accept_prob : 
        the target probability that a step that ends in a new minimum is accepted.
        Use to adjust the temperature
        
        Note: the total acceptance probability is
        
            accrat = 1 - target_new_min_prob * (1 - target_new_min_accept_prob)
        
        so if you want a total acceptance probability of 0.5, you must choose both other 
        probabilities accordingly
        
    interval : 
        the interval at which to adjust temperature and stepsize
    Tfactor : 
        the factor with which to multiply (or divide) the temperature
    sfactor : 
        the factor with which to multiply (or divide) the stepsize
    ediff : 
        if two minima have energies that are within ediff from each other then
        they are considered to be the same minimum
    verbose :
        print status messages
        
    Notes
    -----                    
    We will base the stepsize adjustment on the probability of ending up 
    in a different minimum
    
    We will base the temperature adjustment on the probability of accepting
    a move that ended up in a different minimum
    """

    def __init__(self, stepclass,
                 target_new_min_prob=0.8,
                 target_new_min_accept_prob=0.3,
                 interval=100, Tfactor=0.95, sfactor=0.95, ediff=.001,
                 verbose=False):
        self.stepclass = stepclass
        self.target_new_min_accept_prob = target_new_min_accept_prob
        self.target_new_min_prob = target_new_min_prob
        self.interval = interval
        self.Tfactor = Tfactor
        self.sfactor = sfactor
        self.ediff = ediff
        self.verbose = verbose

        self.energy = None
        self.coords = None

        self.ncalls_tot = 0
        self.reset()

    def reset(self):
        self.nattempts = 0
        self.naccept = 0
        self.nsame = 0


    def takeStep(self, *args, **kwargs):
        """
        basinhopping calls this to take a step
        """
        self.stepclass.takeStep(*args, **kwargs)

    def updateStep(self, accepted, driver=None):
        """
        this is the function basinhopping uses to report results
        """
        self.ncalls_tot += 1
        trial_energy = driver.trial_energy
        if self.energy is None:
            # first time called. Save energy and coords
            self.energy = driver.markovE
            self.coords = driver.coords.copy()
            return

        self.nattempts += 1
        if accepted:
            self.naccept += 1

        # determine if the new minima is the same as the last one
        same = False
        if abs(self.energy - trial_energy) <= self.ediff:
            same = True
            self.nsame += 1

        if not same and accepted:
            self.energy = driver.markovE
            self.coords = driver.coords.copy()

        if self.nattempts % self.interval == 0:
            # if self.verbose:
            #                print "    acceptance probability %.4g" % (float(self.naccept) / self.nattempts)
            self.adjustStep()
            self.adjustTemp(driver)
            self.reset()

    def adjustStep(self):
        """adjust the step size
        
        increase the step size if we're ending up in the same minima too often,
        else decrease the step size
        """
        fnew = 1. - float(self.nsame) / self.nattempts
        if fnew < self.target_new_min_prob:
            self.stepclass.scale(1. / self.sfactor)
        else:
            self.stepclass.scale(self.sfactor)

        # print some status info
        if self.verbose:
            print("adaptive step and temperature: naccept nsame ndiff naccept_diff %d %d %d %d new min probability %.4g" % (
                self.naccept, self.nsame, self.nattempts - self.nsame,
                self.naccept - self.nsame, float(self.naccept - self.nsame) / self.nattempts))
            print("    stepsize    is now %.4g ratio %.4g target %.4g" % (self.stepclass.stepsize,
                                                                          fnew, self.target_new_min_prob))


    def adjustTemp(self, driver):
        """adjust the temperature
        
        increase the temeperature if new minima are rejected too often,
        else decrease the temperature
        """
        ndiff = self.nattempts - self.nsame
        ndiff_accept = self.naccept - self.nsame
        if ndiff == 0:
            faccept = 1
        else:
            faccept = float(ndiff_accept) / ndiff
        if faccept > self.target_new_min_accept_prob:
            driver.acceptTest.temperature *= self.Tfactor
        else:
            driver.acceptTest.temperature /= self.Tfactor
        if self.verbose:
            print("    temperature is now %.4g ratio %.4g target %.4g" % (driver.acceptTest.temperature,
                                                                          faccept, self.target_new_min_accept_prob))


if __name__ == "__main__":
    import numpy as np
    from pele.takestep import displace
    from pele.systems import LJCluster

    natoms = 38
    sys = LJCluster(natoms=38)


    # random initial coordinates
    coords = sys.get_random_configuration()

    takeStep = displace.RandomDisplacement(stepsize=0.4)
    tsAdaptive = AdaptiveStepsizeTemperature(takeStep, interval=300, verbose=True)

    db = sys.create_database()
    opt = sys.get_basinhopping(database=db, takestep=tsAdaptive, coords=coords)
    opt.printfrq = 50
    opt.run(5000)
        

