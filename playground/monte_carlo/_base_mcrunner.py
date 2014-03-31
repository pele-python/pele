import numpy as np
import abc
from pele.optimize import Result
from playground.monte_carlo import MC

class _base_MCrunner(object):
    """
    Abstract method for MC runners, all MC runners should derive from this base class.
    The design of this class relies on a number of implementation choices made for the
    pele::MC cpp class. This is not limiting by any means, developers can easily modify
    this class to write a base class that uses a different basic MC class. 
    Using the pele::MC class is however recommended. 
    *potential should be constructed outside of this class and passed
    *coords are the initial coordinates
    *niter is the total number of MC iterations
    *set_control *MUST* be overwritten in any derived class
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, potential, coords, temperature, stepsize, niter):
        self.potential = potential
        self.ndim = len(coords)
        self.start_coords = coords
        self.temperature = temperature
        self.stepsize = stepsize
        self.niter = niter
        self.mc = MC(self.potential, self.start_coords, self.temperature, self.stepsize)    
        self.result = Result()
        self.result.message = []
    
    @abc.abstractmethod
    def set_control(self, c):
        """set control parameter, this could be temperature or some other control parameter like stiffness of the harmonic potential"""
    
    def get_config(self):
        """Return the coordinates of the current configuration and its associated energy"""
        coords = self.mc.get_coords()
        energy = self.mc.get_energy()
        return coords, energy
    
    def set_config(self, coords, energy):
        """set current configuration and its energy"""
        self.mc.set_coordinates(coords, energy)
    
    def get_results(self):
        """Must return a result object, generally must contain at least final configuration and energy"""
        res = self.result
        res.coords = self.mc.get_coords()
        res.energy = self.mc.get_energy()
        return res
    
    def get_iterations_count(self):
        """returns total number of MCMC steps performed so far"""
        n = self.mc.get_iterations_count()
        return n
    
    def get_accepted_fraction(self):
        """returns the ratio of the number of accepted steps over the total"""
        n = self.mc.get_accepted_fraction()
        return n
    
    def get_stepsize(self):
        """returns the current stepsize"""
        s = self.mc.get_stepsize()
        return s
    
    def get_status(self):
        status = Result()
        status.iteration = self.mc.get_iterations_count()
        status.acc_frac = self.mc.get_accepted_fraction()
        status.conf_reject_frac = self.mc.get_conf_rejection_fraction()
        status.step_size = self.mc.get_stepsize()
        status.energy = self.mc.get_energy()
        status.neval = self.mc.get_neval()
        return status
    
    def run(self):
        """run MCMC walk"""
        self.mc.run(self.niter)