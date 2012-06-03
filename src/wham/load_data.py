import numpy as np 
from numpy import log
#import timeseries # for timeseries analysis 
#import commands
#import pdb;
#import pickle
#import calc_Cv_extrap4_new_utils as whamutil
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *


def exponentialBinEnergy( emin, emax, nebins, Ebar1, sig1, Ebar2, sig2 ):
    """ 
    Calculate bin edges such that dE increase exponentially with E. This is
    so that the lower temperatures that have narrower energy histograms have
    roughly the same number of bins.  If dE increases like

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

    These can be solved analytically for dEmin and dEmax 
    """
    ##########################################################################
    #first determine dEmin and dEmax
    ##########################################################################
    B = (sig2*(Ebar1 - emin) - sig1*(Ebar2-emin))/(sig1-sig2)
    f = (emax - emin)/B + 1.
    dEmin_test = np.log(f) / nebins * B
    dEmax_test = f * dEmin_test
    print "dEmin dEmax ", dEmin_test, dEmax_test
    dEmin = dEmin_test
    dEmax = dEmax_test
    ##########################################################################
    #now generate bin edges
    ##########################################################################
    binenergynp = np.zeros(nebins+1)
    for i in range(nebins):
        binenergynp[i] = emin + dEmin*nebins/log(dEmax/dEmin)*((dEmax/dEmin)**(float(i)/nebins) - 1)
        #print binenergy[i]
    binenergynp[nebins] = emin + dEmin*nebins/log(dEmax/dEmin)*(dEmax/dEmin - 1)
    print "right most bin edge", binenergynp[nebins], "should be ", emax
    return binenergynp


def getNDataLines(fin):
    nlines = 0
    for line in fin:
        if line[0] != '#':
            nlines += 1
    return nlines
    

def readFile(fname, columns):
    data = []
    with open(fname, "r") as fin:
        for line in fin:
            if line[0] == '#': continue
            sline = line.split()
            d = [ float(sline[c]) for c in columns ]
            data.append(d)
    return np.array(data)


class load_data1dExp:
    """
    load data from files and create histograms
    
    filenames: a list of filesnames
    
    ecolumn:   the column in the input files that has the energy
    
    nebins:    number of bins
    
    NEGLECT:   for each replica, discard bins with less fewer visits than 
                NEGLECT*max_visits[i], where max_visits[i] is the maximum number of
                visits for replica i
    
    fskip:     discard the first fraction (fskip) of data
    
    exponential_bins: use exponentially distributed bin sizes. This is so that 
        lower temperatures that have narrower energy histograms have roughly 
        the same number of bins.
    """
    def __init__(self, filenames, ecolumn, nebins = 200, NEGLECT = 0.01, fskip=0., \
                 exponential_bins = True):
    
        nrep = len(filenames)
    
        #loop through files and get maximum and minimum energies
        print "determining emax and emin from input data"
        emin=1e10
        emax=-1e10
        nlines=[]
        Ebar=[]
        Esig=[]
        for f in filenames:
            data=np.genfromtxt(f)
            e=data[:,ecolumn].min()
            if e < emin: emin = e
            e=data[:,ecolumn].max()
            if e > emax: emax = e
            nlines.append(len(data[:,0]))
            ebar = np.mean( data[:,ecolumn] )
            esig = np.std( data[:,ecolumn] )
            Ebar.append(ebar)
            Esig.append(esig)
    
        print "nlines"
        print nlines
    
        emax += 1e-4 #so the highest energy returns the last bin, not one past the last bin
    
        ########################################################
        #determine bin edges
        ########################################################
        self.visits1d = np.zeros([nebins,nrep], np.int32) 
        self.binenergy = np.zeros(nebins)
        binenergynp = np.zeros(nebins+1) #for np.histogram, numpy wants bin endges, including rightmost edge


    
    
        #if exponential_bins then redo the bin calculation
        #use exponentially increasing bin width
        if exponential_bins:
            binenergynp = exponentialBinEnergy( emin, emax, nebins, Ebar[0], Esig[0], Ebar[-1], Esig[-1] )
            self.binenergy[0:nebins] = ( binenergynp[0:-1])
        else:
            dE = (emax-emin)/(nebins)
            print emin, " < E < ", emax, " dE = ", dE
            
            self.binenergy = np.array([ emin + dE*i for i in range(nebins)])
            binenergynp = np.array([ emin + dE*i for i in range(nebins+1)]) #for histogram, numpy wants bin endges, including rightmost edge

        visits1d = self.visits1d

            
    
        #load the energies into visits1d
        print "reading data from files"
        for k in range(nrep):
            with open(filenames[k],"r") as f:
                lcount = -1
                lfirst = nlines[k]*fskip
                print "lfirst ", lfirst
                edata=np.genfromtxt(f)[lfirst:,ecolumn]
                visits1d[:,k], retbins = np.histogram(edata, binenergynp)
    
        #set visits to zero if the count is below a threshold
        #the threshold is different for each replica, and is determined by the ratio to the maximum count
        for k in range(nrep):
            visitsmax=visits1d[:,k].max()
            ind=np.where(visits1d.astype(float)[:,k]/visitsmax<NEGLECT)
            visits1d[ind[0],k] = 0

def determineBinEdge(nbins, datalist, column, minmax = [], exponential_bins=True):
    """
    return the bin edges, including the rightmost edge
    """
    #get min and max
    if len(minmax) == 2:
        emin = minmax[0]
        emax = minmax[1]
    else:
        emin = min( np.min( data[:,column] ) for data in datalist )
        emax = max( np.max( data[:,column] ) for data in datalist )
    
    print "min max", emin, emax
    emax += 1e-4 #so the highest energy returns the last bin, not one past the last bin

    if exponential_bins:
        Ebar = [ np.mean(data[:,column]) for data in datalist ]
        Esig = [ np.std(data[:,column]) for data in datalist ]
        binenergynp = exponentialBinEnergy( emin, emax, nbins, Ebar[0], Esig[0], Ebar[-1], Esig[-1] )
    else:
        dE = (emax-emin)/(nbins)
        binenergynp = np.array([ emin + dE*i for i in range(nbins+1)]) 
        print emin, " < E < ", emax, " dE = ", dE
    
    return binenergynp
        



def loadData2d(filenames, columns = [0], fskip=0., qcombine=[]):
    """
    load data from columns in filenames.  return the data as a list of numpy arrays
    
    fskip:  discard the first fraction (fskip) of lines
    
    qcombine:  optional:  qcombine = [qcolumn1, qcolumn2, frac]
        make a linear combination of order parameters from columns 1 and 2
        q = (q1 + q2weight*q2)/(1+q2weight)
        
        if this is chosen, the columns will be appended to the other given columns

    """
    #qcombine = [qcolumn1, qcolumn2, frac]
    #make a linear combination of order parameters from columns 1 and 2
    # q = (q1 + q2weight*q2)/(1+q2weight)
    qcomb = False
    if len(qcombine) == 3:
        qcomb = True
        qcolumn1 = int(qcombine[0])
        qcolumn2 = int(qcombine[1])
        q2weight = qcombine[2]   
        print "qcombine: ", qcolumn1, qcolumn2, q2weight
        columns.append( qcolumn1 ) 
        columns.append( qcolumn2 )

    #if the columns have a header, print the name of the column requested
    with open(filenames[0], "r") as fin:
        line = fin.readline()
        sline = line.split()
        if sline[0][0] == '#':
            print "column headers ",
            for c in columns:
                print sline[c],
            print


    #######################################################
    #load data
    #######################################################
    #load the data
    print "loading data"
    mydatalist = []
    for fname in filenames:
        mydata = readFile(fname, columns)
        mydatalist.append(mydata)
        print "mydata", np.shape(mydata)

    #discard the first fraction (fskip) lines
    if fskip > 0:
        print "fskip", fskip
        mydatalistnew = []
        for data in mydatalist:
            imin = int(fskip * len(data[:,0]))
            #print "first index", imin
            #data = data[imin:,:]
            datanew = data[imin:,:]
            mydatalistnew.append(datanew)
        mydatalist = mydatalistnew


    if qcomb:
        #recombine data
        mydatanew = []
        for data in mydata:
            datanew = np.zeros([ len(data[:,0]), 2] )
            datanew[:,0] = data[:,0]
            datanew[:,0] = (data[:,1] + data[:,2]*q2weight ) / (1.+q2weight)
            mydatanew.append( data  )
        mydatalist = mydatanew
    return mydatalist


def binData2d( binenergynp, binqnp, datalist, NEGLECT = 0.01):

    nebins = len(binenergynp) -1
    nqbins = len(binqnp) - 1
    nreps = len(datalist)
    
    visits2d = np.zeros([nebins,nqbins,nreps], np.integer) 
        

    #load the energies into visits2d
    for k, data in enumerate(datalist):
        e = data[:,0]
        q = data[:,1]
        visits2d[:,:,k], xbins, ybins = np.histogram2d(e, q, [binenergynp, binqnp])


    #set visits to zero if the count is below a threshold
    #the threshold is different for each replica, and is determined by the ratio to the maximum count
    for k in range(nreps):
        visitsmax=visits2d[:,:,k].max()
        minvis = max(1, int(visitsmax * NEGLECT))
        ind = np.where(visits2d[:,:,k] <= minvis)
        print k, "total visits ", visits2d[:,:,k].sum(), "max visits ", visitsmax, "max visits to be discarded ", visits2d[ind[0],ind[1],k].max()
        visits2d[ind[0],ind[1],k] = 0
    
    return visits2d

