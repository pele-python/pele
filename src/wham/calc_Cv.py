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
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.pyplot import *
#mbar = pickle.load(open("mbar.pickle","rb"))

#parser = argparse.ArgumentParser(description='Combine energy and overlap data from multiple runs at different temperatures into one histogram and print F_q.')
#parser.add_argument('-f', type=int, nargs='1', help='the number of free particles')
#parser.parse_args()

def usage():
    print sys.argv[0], " [-hF -o output_prefix -r rskip -q qcolumn -e ecolumn -E nebins -c input -T TRANGE] -f nfree"
    print 'Combine energy and from multiple runs at different temperatures into one histogram and print various quantities.'
    print '  -h print this help and exit'
    print '  -f nfree : number of mobile particles, used to determine nqbins and # degrees of freedom'
    print '  -o output_prefix : change the default output_prefix'
    print '  -F :  dont use pickle file'
    print '  -r rskip : skip the first fraction r of the data files'
    print '  -q qcolumn : which column to get the overlap data from'
    print '  -e ecolumn : which column to get the energy data from'
    print '  -E nebins  : number of energy bins (=300)'
    print '  -c input : Make a linear combination of two order parameters.'
    print '             Input will have the form "q1column q2column q2weight"'
    print '             The order parameter will be q = (q1 + q2weight*q2)/(1+q2weight)'
    print '  -T TRANGE : set TRANGE for the calculation of Fq.  TRANGE should have the format "Tmin Tmax numT"'

usepkl=True
nfree=1;
outprefix="out"
rskip=0.0
qcolumn=3
ecolumn=2
qcombine=[]
TRANGEi=[]
nebins=300


opts, args = getopt.getopt(sys.argv[1:], "hf:o:Fr:q:c:T:e:E:", ["help", "nfree="])
output = None
verbose = False
for o, a in opts:
    if o == "-f":
        nfree=int(a)
    elif o == "-o":
        outprefix=a
        print "output_prefix = ", outprefix
    elif o == "-F":
        usepkl=False
    elif o == "-r":
        rskip=float(a)
        print "will skip the first ", rskip, "of the data files"
    elif o == "-q":
        qcolumn=int(a)
        print "using qcolumn = ", qcolumn
    elif o == "-e":
        ecolumn=int(a)
        print "using ecolumn = ", ecolumn
    elif o == "-E":
        nebins=int(a)
        print "using nebins = ", nebins
    elif o == "-c":
        qcombline=a
        qcombine = [float(b) for b in qcombline.split()]
        if len(qcombine) != 3:
            print "-c: qcombine must have 3 parts: ", qcombine
            usage()
            exit(1)
        print "using qcombine: ", qcombine
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
        assert False, "unhandled option"



pklname=outprefix+".pickle"
if not usepkl or not os.path.isfile(pklname):

    #get temperatures and filenames
    Tlist=list(np.genfromtxt('temperatures'))
    rep1=0 #don't use the first 8 replicas.  There is a better way to do this.
    Tlist=Tlist[(rep1):]
    nrep = len(Tlist)
    filenames=['overlap.'+str(n+rep1+1) for n in range(nrep)]
    #ecolumn=2
    #qcolumn=3

    #OK, now we have a list of temperatures and filenames for each replicas
    print "replica list:"
    for n in range(nrep):
        print Tlist[n], filenames[n]


    #nfree=64; #THIS SHOULD BE PASSED TO THE PROGRAM
    print "USING nfree = ", nfree
    #nebins = 300
    nqbins=nfree+1
    #nrep=len(Tlist)
    #qmin=0.
    #qmax=1.


    data = load_data.load_data1dExp(filenames, ecolumn, nebins=nebins, fskip=rskip)
    visits1d = data.visits1d

    wham = WHAM.wham1d(Tlist, data.binenergy, visits1d)

    wham.minimizeNew()
    #wham.minimize()

    print "dumping WHAM1d to pickle file: ", pklname
    pickle.dump(wham,open(pklname,"wb"))

else:
    print "=================================================================="
    print "loading WHAM2d from pickle file: ", pklname
    print "=================================================================="
    wham = pickle.load(open(pklname,"rb")) 




fname=outprefix+".Cv"
NDOF = nfree*3
print 'writing Cv to', fname, ': ', NDOF, "degrees of freedom"
fout=open(fname,"w")
wham.calc_Cv(NDOF, fout)
fout.close()

fname=outprefix+".n_E"
print 'writing density of states, log(n_E) to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if wham.visits1d[i,:].sum() > 0:
            val=wham.logn_E[i]
            fout.write( str(wham.binenergy[i]) +" "+ str(val) + "\n" ) 

fname=outprefix+".visits"
print 'writing sum visits to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if wham.visits1d[i,:].sum() > 0:
            val = wham.visits1d[i,:].sum()
            fout.write( str(wham.binenergy[i]) +" "+ str(val) + "\n" ) 

fname=outprefix+".visits.all"
print 'writing visits to', fname
with open(fname,"w") as fout:
    for i in range((wham.nebins)):
        if wham.visits1d[i,:].sum() > 0:
            vallist = (wham.visits1d[i,:].tolist())
            slist = [ str(v) for v in vallist]
            fout.write( str(wham.binenergy[i]) +" ".join(slist) + "\n" ) 
