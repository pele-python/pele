import numpy as np
from numpy import log, exp



def logSum1(a, b):
    """
    return log( exp(a) + exp(b) )
    """
    if a > b:
        return a + log(1.0 + exp(-a + b) )
    else:
        return b + log(1.0 + exp(a - b) )


def logSum(log_terms):
    """
    Compute the log of a sum of terms whose logarithms are provided.

    REQUIRED ARGUMENTS  
      log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

    RETURN VALUES
      log_sum is the log of the sum of the terms.
    """
    try:
        return logSumFast(log_terms)
    except ImportError:
        #print "using slow logsum, install scipy weave inline blitz to speed up"
        return logSumSlow(log_terms)

def logSumSlow(log_terms):
    log_sum = log_terms[0] 
    for lt in log_terms[1:]:
        if log_sum > lt:
            log_sum = log_sum + log(1.0 + exp(-log_sum + lt) )
        else:
            log_sum = lt + log(1.0 + exp(log_sum - lt) )
    # return the log sum
    return log_sum

def logSumFort(log_terms):
    import _wham_utils
    return _wham_utils.logsum(log_terms)

def logSumWeave(log_terms):
    from scipy import weave
    from scipy.weave import converters
    nterms = len(log_terms)
    code = """
    double log_sum = log_terms(0);
    for (int i = 1; i < nterms; ++i){
        double lt = log_terms(i);
        if (log_sum > lt){
            log_sum = log_sum + log(1.0 + exp(-log_sum + lt) );
        } else {
            log_sum = lt + log(1.0 + exp(log_sum - lt) );
        }
    }
    return_val = log_sum;
    """
    log_sum = weave.inline(code, ["nterms", "log_terms"], type_converters=converters.blitz, verbose=2)
    return log_sum

def logSumFast(log_terms):
    """
    there are two options available in case the user doesn't have package scipy (as
    is the case on some clusters) or doesn't have f2py
    """
    try:
        return logSumFort(log_terms)
    except ImportError:
        return logSumWeave(log_terms)

def calc_Cv(logn_E, visits1d, binenergy, NDOF, Treplica, k_B, TRANGE=None, NTEMP=100, use_log_sum = None):
    if use_log_sum == None:
        try:
            logSumFast( np.array([.1, .2, .3]) )
            use_log_sum = True
            #print "using logsum"
        except ImportError:
            #dont use log sum unless the fast logsum is working
            #print "not using logsum because it's too slow.  Install scipy weave to use the fast version of logsum"
            use_log_sum = False
    else:
        #print "use_log_sum = ", use_log_sum
        pass
    

    #put some variables in this namespace
    nrep, nebins = np.shape(visits1d)
    #print "nreps, nebins", nrep, nebins

    nz = np.where( visits1d.sum(0) != 0)[0]

    #find the ocupied bin with the minimum energy
    EREF = np.min(binenergy[nz]) - 1.

    #now calculate partition functions, energy expectation values, and Cv
    if TRANGE == None:
        #NTEMP = 100 # number of temperatures to calculate expectation values
        TMAX = Treplica[-1]
        TMIN = Treplica[0]
        TINT=(TMAX-TMIN)/(NTEMP-1)
        TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]

    dataout = np.zeros( [NTEMP, 6] )
    if abs((binenergy[-1] - binenergy[-2]) - (binenergy[-2] - binenergy[-3]) ) > 1e-7:
        print "calc_Cv: WARNING: dE is not treated correctly for exponential energy bins"

    for count,T in enumerate(TRANGE):
        kBT = k_B*T
        #find expoffset so the exponentials don't blow up
        dummy = logn_E[nz] - (binenergy[nz] - EREF)/kBT
        expoffset = np.max(dummy)
        if not use_log_sum:
            dummy = np.exp(dummy)
            Z0 = np.sum(dummy)
            Z1 = np.sum( dummy *  (binenergy[nz] - EREF) )
            Z2 = np.sum( dummy *  (binenergy[nz] - EREF)**2 )
            lZ0 = np.log(Z0)
            lZ1 = np.log(Z1)
            lZ2 = np.log(Z2)
        else:
            lZ0 = logSum( dummy )
            lZ1 = logSum( dummy + log(binenergy[nz] - EREF) )
            lZ2 = logSum( dummy + 2.*log(binenergy[nz] - EREF) )
            Z0 = np.exp( lZ0 )
            Z1 = np.exp( lZ1 )
            Z2 = np.exp( lZ2 )

        i = nebins-1
        #if abs(binenergy[-1] - binenergy[-2]) > 1e-7:
        #    print "calc_Cv: WARNING: dE is not treated correctly for exponential bins"
        if i == nebins-1:
            dE = binenergy[i]-binenergy[i-1]
        else:
            dE = binenergy[i+1]-binenergy[i]
            
        if dE/kBT < 1.0E-7:
            ONEMEXP=-dE/kBT
        else:
            ONEMEXP= 1.0-np.exp(dE/kBT)

        Eavg = NDOF*kBT/2.0 + 1.0*(kBT + dE/ONEMEXP) + exp(lZ1-lZ0) + EREF
        
        Cv = NDOF/2. + 1.0*(1.0 - dE**2 * exp(dE/kBT)/(ONEMEXP**2*kBT**2)) \
            - exp(lZ1 - lZ0)**2 / kBT**2 + exp(lZ2-lZ0) / kBT**2
            #- (Z1/(Z0*kBT))**2 + Z2/(Z0*kBT**2)
        
        dataout[count,0] = T
        dataout[count,1] = lZ0 + expoffset
        dataout[count,2] = lZ1 + expoffset
        dataout[count,3] = lZ2 + expoffset
        dataout[count,4] = Eavg
        dataout[count,5] = Cv

        
        #np.array([kBT, Z0*np.exp(expoffset), Z1*np.exp(expoffset), Z2*np.exp(expoffset), Eavg, Cv, np.log(Z0)+expoffset, np.log(Z1)+expoffset, np.log(Z2)+expoffset]).tofile(fout," ")
        
        #fout.write("\n")

    return dataout


import unittest
class TestLogSum(unittest.TestCase):
    def testLogSum(self):
        vals = np.random.rand(50) + 0.001
        ls1 = logSumSlow(vals)
        ls2 = logSumWeave(vals)
        print "%g - %g = %g" % (ls1, ls2, ls1-ls2)
        self.assertTrue( abs(ls1 - ls2) < 1e-12, "logSumFast is different from logSumSlow: %g - %g = %g" % (ls1, ls2, ls1-ls2) )
    def testLogFort(self):
        vals = np.random.rand(50) + 0.001
        ls1 = logSumSlow(vals)
        ls2 = logSumFort(vals)
        print "%g - %g = %g" % (ls1, ls2, ls1-ls2)
        self.assertTrue( abs(ls1 - ls2) < 1e-12, "logSumFast is different from logSumSlow: %g - %g = %g" % (ls1, ls2, ls1-ls2) )

        

if __name__ == "__main__":
    unittest.main()
