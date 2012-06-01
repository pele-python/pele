from numpy import *
import numpy as np #to access np.exp() not built int exp
from math import *
import scipy.optimize
#import timeseries # for timeseries analysis 
import commands
import pdb;
import pickle
#import calc_Cv_extrap4_new_utils as whamutil
from wham_potential import WhamPotential
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *


def exponentialBinEnergy( binenergynp, emin, emax, nebins, Ebar1, sig1, Ebar2, sig2 ):
    """ 
    Caluclate bin edges such that dE increase exponentially with E. This is
    so that the lower temperatures that have narrower energy histograms have
    roughly the same number of bins.  If dE increses like

        dE[i] = dEmin * (dEmax/dEmin)**(i/nebins)

    then the energy of bin i is

        binenergy[i] = emin + dEmin*nebins / log(dEmax/dEmin) * (exp( i/nebins * log(dEmax/dEmin) ) - 1)

    We have two free parameters, dEmax and dEmin which are determined as follows.  Let

        sig1 = std. dev. of the energy of the lowest temp
        Ebar1 = mean of the energy of the lowest temp
        sig2 = std. dev. of the energy of the highest temp
        Ebar2 = mean of the energy of the highest temp

    Then the two conditions for determining dEmin and dEmax are

        dE( i(Ebar1) ) / sig1  = dE( i(Ebar2) ) / sig2
        E( nebins ) = emax

    These can be solved analitically for dEmin and dEmax 
    """
    ##########################################################################
    #first determine dEmin and dEmax
    ##########################################################################
    B = (sig2*(Ebar1 - emin) - sig1*(Ebar2-emin))/(sig1-sig2)
    f = (emax - emin)/B + 1.
    dEmin_test = log(f) / nebins * B
    dEmax_test = f * dEmin_test
    print "dEmin dEmax ", dEmin_test, dEmax_test
    dEmin = dEmin_test
    dEmax = dEmax_test
    ##########################################################################
    #now generate bin edges
    ##########################################################################
    for i in range(nebins):
        binenergynp[i] = emin + dEmin*nebins/log(dEmax/dEmin)*((dEmax/dEmin)**(float(i)/nebins) - 1)
        #print binenergy[i]
    binenergynp[nebins] = emin + dEmin*nebins/log(dEmax/dEmin)*(dEmax/dEmin - 1)
    print "right most bin edge", binenergynp[nebins], "should be ", emax



class load_data1dExp:
  def __init__(self, filenames, ecolumn, nebins = 200, NEGLECT = 0.01, fskip=0.):

    nrep = len(filenames)


    #loop through files and get maximum and minimum energies
    print "determining emax and emin from input data"
    emin=1e10
    emax=-1e10
    nlines=[]
    Ebar=[]
    Esig=[]
    for f in filenames:
        data=genfromtxt(f)
        e=data[:,ecolumn].min()
        if e < emin: emin = e
        e=data[:,ecolumn].max()
        if e > emax: emax = e
        nlines.append(len(data[:,0]))
        ebar = mean( data[:,ecolumn] )
        esig = std( data[:,ecolumn] )
        Ebar.append(ebar)
        Esig.append(esig)

    print "nlines"
    print nlines

    emax += 1e-4 #so the highest energy returns the last bin, not one past the last bin

    dE = (emax-emin)/(nebins)
    print emin, " < E < ", emax, " dE = ", dE

    self.visits1d = zeros([nebins,nrep], int32) 
    self.binenergy = array([ emin + dE*i for i in range(nebins)])
    visits1d = self.visits1d
    binenergy = self.binenergy
    binenergynp = array([ emin + dE*i for i in range(nebins+1)]) #for histogram, numpy wants bin endges, including rightmost edge

    #use exponentially increasing bin width
    exponential_bins = True
    if exponential_bins:
        exponentialBinEnergy( binenergynp, emin, emax, nebins, Ebar[0], Esig[0], Ebar[-1], Esig[-1] )
        binenergy[0:nebins] = ( binenergynp[0:-1])

    #load the energies into visits1d
    print "reading data from files"
    for k in range(nrep):
        with open(filenames[k],"r") as f:
            lcount = -1
            lfirst = nlines[k]*fskip
            print "lfirst ", lfirst
            edata=genfromtxt(f)[lfirst:,ecolumn]
            visits1d[:,k], retbins = np.histogram(edata, binenergynp)

    #set visits to zero if the count is below a threshold
    #the threshold is different for each replica, and is determined by the ratio to the maximum count
    for k in range(nrep):
        visitsmax=visits1d[:,k].max()
        ind=where(visits1d.astype(float)[:,k]/visitsmax<NEGLECT)
        visits1d[ind[0],k] = 0




