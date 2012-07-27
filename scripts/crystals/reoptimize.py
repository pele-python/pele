import pickle
import dmagmin_ as GMIN
from pygmin.utils import crystals
from pygmin.potentials import gminpotential 
from optparse import OptionParser

# add some program options
parser = OptionParser(usage = "usage: %prog [options] arg1 arg2")
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("--tol",
                  dest="tol", default=1e-4, action="store",type="float",
                  help="tolerance for quench")
parser.add_option("--maxErise",
                  dest="maxErise", default=1e-2, action="store",type="float",
                  help="maximum increase in energy for lbfgs")
parser.add_option("-o","--out",
                  dest="out",# action="store",type="float",
                  help="store in this new object. If not given, the input file is overwritten")

(options, args) = parser.parse_args()

# print help if no input file is given
if(len(args) != 1):
    parser.print_help()
    exit(-1)

infile = args[0]
outfile = infile
# use input file for output unless specified otherwise
if(options.out!=None):
    outfile = options.out
   
# fake compareMinima function
def compareMinima():
    pass

# initialize GMIN
GMIN.initialize()
pot = gminpotential.GMINPotental(GMIN)
crystals.GMIN = GMIN

# open the storage class
save = pickle.load(open(infile,"r"))

#loop over minima
i=0
for m in save.data:
    i+=1
    ret = crystals.quenchCrystal(m.coords, pot.getEnergyGradient, tol=options.tol, maxErise=options.maxErise)
    print "minimum",i,"energy (old energy)",ret[1],"(%f)"%(m.E)
    m.E = ret[1]
    m.coords = ret[0]
    pickle.dump(save, open(outfile, "w"))
    
