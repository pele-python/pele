import dmagmin_ as GMIN
import pickle
#import sys

from optparse import OptionParser

# add some program options
parser = OptionParser(usage = "usage: %prog [options] storage")
parser.add_option("--cif",
                  dest="writeCIF", action="store_true",
                  help="export cif files")
parser.add_option("--cif-dir",
                  dest="cifDir", default=".", action="store",type="string",
                  help="directory to write cifs to")

(options, args) = parser.parse_args()

# print help if no input file is given
if(len(args) != 1):
    parser.print_help()
    exit(-1)

def compareMinima():
    pass

if(options.writeCIF):
    GMIN.initialize()

save = pickle.load(open(args[0]))

i=0
for m in save.data:
    i+=1
    print "minimum",i, "energy",m.E
    if(options.writeCIF):
        filename = options.cifDir+"/lowest%03d.cif"%(i)
        print "writing",filename
        GMIN.writeCIF(filename, m.coords, "E"+str(m.E))
