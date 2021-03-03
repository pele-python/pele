"""
Create and print histograms.  Especially energy histograms.

.. currentmodule:: pele.utils.histogram

.. autosummary::
    :toctree: generated
    
    EnergyHistogram
    PrintHistogram
"""
from __future__ import print_function

import numpy as np

__all__ = ["EnergyHistogram", "PrintHistogram"]


class EnergyHistogram(object):
    """
    this class will build 1 dimensional histogram.  
    It's designed for energies, but it could work for any float data
    """

    def __init__(self, emin, emax, nbins=100):
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        self.de = (self.emax - self.emin) / self.nbins

        self.visits = np.zeros(self.nbins)
        self.count = 0

    def insert(self, e):
        if not self.emin <= e < self.emax:
            print("histogram> warning: energy out of range", e)
            return
        i = int((e - self.emin) / self.de)
        self.visits[i] += 1
        self.count += 1

    def insertWrapper(self, e, coords, acceptstep):
        return self.insert(e)

    def __iter__(self):
        return HistIter(self)

    def __call__(self, e, coords, acceptstep):
        return self.insert(e)


class HistIter(object):
    def __init__(self, hist):
        self.hist = hist
        self.counter = -1

    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        if self.counter >= self.hist.nbins:
            raise StopIteration
        return (self.hist.emin + self.hist.de * self.counter), \
               self.hist.visits[self.counter]

    def next(self):
        return self.__next__()


class PrintHistogram(object):
    def __init__(self, fname, hist, interval):
        self.fname = fname
        self.hist = hist
        self.interval = interval
        self.outstream = open(self.fname, "w")
        self.count = 0

    def printHis(self):
        for e, count in self.hist:
            # print "%g %d\n" % (e, count),
            self.outstream.write("%g %d\n" % (e, count))
        self.outstream.write("\n\n")


    def printEvent(self):
        self.count += 1
        if self.count % self.interval == 0:
            self.printHis()

    def __call__(self, a, b, c, **kwargs):
        self.printEvent()
            

