#!/usr/bin/python
import numpy as np #to access np.exp() not built int exp
#import timeseries # for timeseries analysis 
#import commands
#import pdb;
import pickle
import os.path
#import calc_Cv_extrap4_new_utils as whamutil
import histogram_reweighting1d as WHAM
import load_data
#import utils
import getopt, sys
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.pyplot import *
#mbar = pickle.load(open("mbar.pickle","rb"))

#parser = argparse.ArgumentParser(description='Combine energy and overlap data from multiple runs at different temperatures into one histogram and print F_q.')
#parser.add_argument('-f', type=int, nargs='1', help='the number of free particles')
#parser.parse_args()

def usage():
    print sys.argv[0], " [-hF -o output_prefix -r rskip -q qcolumn -e ecolumn -E nebins -c input -T TRANGE] -f nfree"
    print 'Combine energy from multiple runs at different temperatures into one histogram and print various quantities.'
    print 'The temperatures of the replicas should be stored in file "temperatures"'
    print '  -h               : print this help and exit'
    print '  -f nfree         : number of mobile particles used to determine number of degrees of freedom'
    print '  -o output_prefix : change the default output_prefix'
    print '  -i input_prefix  : data should be in files of form input_prefix.#'
    print '  -F               :  dont use pickle file'
    print '  -r rskip         : skip the first fraction r of the data files'
    print '  -e ecolumn       : which column to get the energy data from'
    print '  -E nebins        : number of energy bins (=300)'
    print '  -T TRANGE        : set TRANGE for the calculation of Cv.  '
    print '                     TRANGE should have the format "Tmin Tmax numT" (not implemented)'

usepkl=True
nfree=1;
outprefix="out"
inprefix = "overlap"
rskip=0.0
ecolumn=2
qcombine=[]
TRANGEi=[]
nebins=300



opts, args = getopt.getopt(sys.argv[1:], "hf:o:Fr:T:e:E:i:", ["help", "nfree="])
output = None
verbose = False
for o, a in opts:
    if o == "-f":
        nfree=int(a)
    elif o == "-o":
        outprefix=a
        print "output_prefix = ", outprefix
    elif o == "-i":
        inprefix=a
        print "output_prefix = ", inprefix
    elif o == "-F":
        usepkl=False
    elif o == "-r":
        rskip=float(a)
        print "will skip the first ", rskip, "of the data files"
    elif o == "-e":
        ecolumn=int(a)
        print "using ecolumn = ", ecolumn
    elif o == "-E":
        nebins=int(a)
        print "using nebins = ", nebins
    elif o == "-T":
        line = [float(b) for b in a.split()]
        if len(line) != 3:
            print "-T: TRANGE must have 3 parts: ", qcombine
            usage()
            exit(1)
        TMIN = line[0]
        TMAX = line[1]
        NTEMP = int(line[2])
        TINT=(TMAX-TMIN)/(NTEMP-1)
        TRANGEi = [ TMIN + i*TINT for i in range(NTEMP) ]
        print "using TRANGE: ", TRANGEi
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    else:
        print "don't understand argument", o, a
        usage()
        assert False, "unhandled option"



pklname=outprefix+".pickle"
if not usepkl or not os.path.isfile(pklname):

    #get temperatures and filenames
    Tlist=list(np.genfromtxt('temperatures'))
    rep1=0 #don't use the first rep1 replicas.  There is a better way to do this.
    Tlist=Tlist[(rep1):]
    nrep = len(Tlist)
    filenames=[inprefix+'.'+str(n+rep1+1) for n in range(nrep)]

    #OK, now we have a list of temperatures and filenames for each replicas
    print "replica list:"
    for n in range(nrep):
        print Tlist[n], filenames[n]


    print "USING nfree = ", nfree


    #load data
    datalist = load_data.loadData(filenames, [ecolumn], fskip=rskip )
    
    #determine bin edges
    binenergy1 = load_data.determineBinEdge(nebins, datalist, column=0, exponential_bins=True)
    
    #histogram the data
    visits1d = load_data.binData1d(binenergy1, datalist)
    #visits1d = np.transpose(visits1d)

    wham = WHAM.wham1d(Tlist, binenergy1[:-1], visits1d)

    wham.minimize()
    #wham.globalMinimization()

    print "dumping WHAM1d to pickle file: ", pklname
    pickle.dump(wham,open(pklname,"wb"))

else:
    print "=================================================================="
    print "loading WHAM from pickle file: ", pklname
    print "=================================================================="
    wham = pickle.load(open(pklname,"rb")) 


visits1d = np.transpose(wham.visits1d)


fname=outprefix+".Cv"
NDOF = nfree*3
print 'writing Cv to', fname, ': ', NDOF, "degrees of freedom"
cvdata = wham.calc_Cv(NDOF)
with open(fname,"w") as fout:
    np.savetxt(fout, cvdata )
fname += ".pdf"
print 'saving Cv plot to', fname
with open(fname,"w") as fout:
    plt.xlabel("T")
    plt.ylabel("Cv")
    plt.plot(cvdata[:,0], cvdata[:,5], "-")
    plt.savefig(fout)


    

fname=outprefix+".n_E"
print 'writing density of states, log(n_E) to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if visits1d[i,:].sum() > 0:
            val=wham.logn_E[i]
            fout.write( str(wham.binenergy[i]) +" "+ str(val) + "\n" ) 

fname=outprefix+".visits"
print 'writing sum visits to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if visits1d[i,:].sum() > 0:
            val = visits1d[i,:].sum()
            fout.write( str(wham.binenergy[i]) +" "+ str(val) + "\n" ) 

fname=outprefix+".visits.all"
print 'writing visits to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if visits1d[i,:].sum() > 0:
            vallist = (visits1d[i,:].tolist())
            slist = [ str(v) for v in vallist]
            fout.write( str(wham.binenergy[i]) +" ".join(slist) + "\n" ) 
