import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic
from playground.monte_carlo import RandomCoordsDisplacement, MetropolisTest, RecordEnergyHistogram, MC

def test_mc():
    niter = 10000000
    ndim = 5
    temperature = 1.0
    stepsize = 0.8/np.sqrt(ndim)
    bin = 0.01
    k = 1.0
    Emax = 5.0
    Emin = 0.0
    
    #declare np array of desired dimensionality for origin of harmonic potential
    origin = np.zeros(ndim)
    #declare np array for actual coordinates
    #coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #must put coordinates sampled from Pow(ndim)
    coords = np.zeros(ndim)
    #declare np array where to save histogram at the end
    
    harmonic = Harmonic(origin,k)
    histogram = RecordEnergyHistogram(0,30,bin)
    step = RandomCoordsDisplacement(ndim)
    metropolis = MetropolisTest()
    mc = MC(harmonic, coords, temperature, stepsize)
    mc.set_takestep(step)
    mc.add_accept_test(metropolis)
    mc.add_action(histogram)
    
    mc.run(niter)
    
    #histogram.print_histogram(niter)
    n = histogram.get_histogram_size()
    hist = np.zeros(n)
    histogram.get_histogram(hist)
    
#    for i in xrange(len(hist)):
#        print hist[i]
        
    # Open a file
    myfile = open("harmonic_dos.dat", "wb")
    
    E=0
    sumh = 0

    for i in xrange(len(hist)):
        if E >= Emin and E <= Emax:
            sumh += hist[i] * np.exp(temperature*E) * bin;
            E += bin;
          
    E = 0
    
    for i in xrange(len(hist)):
        h = hist[i] * np.exp(temperature*E) / sumh;
        myfile.write("{0} \t {1} \n".format(E,h))
        E += bin;
    
    myfile.close()
    
    myfile = open("analytical_dos.dat", "wb")
    
    E=0
    sumh = 0

    for i in xrange(len(hist)):
        if E >= Emin and E <= Emax:
            sumh += pow(E,ndim/2-1) * bin;
            E += bin;
          
    E = 0
    
    for i in xrange(len(hist)):
        h = pow(E,ndim/2-1)/sumh;
        myfile.write("{0} \t {1} \n".format(E,h))
        E += bin;
    
    myfile.close()

if __name__ == "__main__":
    test_mc()

