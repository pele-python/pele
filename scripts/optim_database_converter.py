import argparse
import os
import numpy as np

from pele.storage import Database
from pele.utils.optim_compatibility import OptimDBConverter    

    

def main():
    parser = argparse.ArgumentParser(description="""
convert an OPTIM database to a pele sqlite database.  Four files are needed.  Normally they are called:

    points.min : the coordinates of the minima in binary format
    min.data   : additional information about the minima (like the energy)
    points.ts  : the coordinates of the transition states
    min.ts     : additional information about transition states (like which minima they connect)

Other file names can optionally be passed.  Some fortran compilers use non-standard endianness to save the
binary data.  If your coordinates are garbage, try changing the endianness.
    """, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--ndof', help='Number of total degrees of freedom (e.g. 3*number of atoms).  This is simply the length of a coordinates vector.', 
                        type=int, default=None)
    parser.add_argument('--Database','-d', help = 'Name of database to write into', type = str, default="optimdb.sqlite")
    parser.add_argument('--Mindata','-m', help = 'Name of min.data file', type = str, default="min.data")
    parser.add_argument('--Tsdata','-t', help = 'Name of ts.data file', type = str, default="ts.data")
    parser.add_argument('--Pointsmin','-p', help = 'Name of points.min file', type = str, default="points.min")
    parser.add_argument('--Pointsts','-q', help = 'Name of points.ts file', type = str, default="points.ts")
    parser.add_argument('--endianness', help = 'set the endianness of the binary data.  Can be "<" for little-endian or ">" for big-endian', type = str, default="=")
    args = parser.parse_args()
    
    db = Database(args.Database)
    
    cv = OptimDBConverter(database=db, ndof=args.ndof, mindata=args.Mindata, 
                 tsdata=args.Tsdata, pointsmin=args.Pointsmin, pointsts=args.Pointsts,
                 endianness=args.endianness)

    cv.setAccuracy()
     
    cv.convert()
    cv.db.session.commit()


    
if __name__ == "__main__":
    main()

