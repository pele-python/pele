import pickle
import dmagmin_ as GMIN
from pele.utils import crystals,dmagmin
from pele.potentials import gminpotential 
from optparse import OptionParser
from pele.storage.database import Database

# add some program options
parser = OptionParser(usage = "usage: %prog [options] database")
#parser.add_option("-f", "--file", dest="filename",
#                 help="write report to FILE", metavar="FILE")
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
   
# initialize GMIN
GMIN.initialize()
pot = gminpotential.GMINPotential(GMIN)
crystals.GMIN = GMIN

# open the storage class
db_in = Database(db=infile)
if(options.out!=None):
    outfile = options.out
    db_out = Database(db=outfile)

print "Start to reoptimize minima"
for m in db_in.minima():
    print "optimizing minima",m._id,", energy", m.energy
    ret = dmagmin.quenchCrystal(m.coords, pot.getEnergyGradient, tol=options.tol, maxErise=options.maxErise)
    print "new energy",ret[1],"(difference %f)"%(ret[1] - m.energy)
    print
    if(options.out!=None):
        db_out.addMinimum(ret[1], ret[0])
    else:
        m.energy = ret[1]
        m.coords = ret[0]
        db_in.session.commit()
    
