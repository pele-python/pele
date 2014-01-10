import argparse
import numpy as np

from pele.gui import run_gui
from pele.storage import Database
from pele.systems import BaseSystem

from optim_database_converter import OptimDBConverter

        

def main():
    parser = argparse.ArgumentParser(description="start a limited functionality gui from an OPTIM database", 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('ndof', help='Number of total degrees of freedom (e.g. 3*number of atoms).  This is simply the length of a coordinates vector.', 
                        type=int)
#    parser.add_argument('--Database','-d', help = 'Name of database to write into', type = str, default="optimdb.sqlite")
    parser.add_argument('--Mindata','-m', help = 'Name of min.data file', type = str, default="min.data")
    parser.add_argument('--Tsdata','-t', help = 'Name of ts.data file', type = str, default="ts.data")
    parser.add_argument('--Pointsmin','-p', help = 'Name of points.min file', type = str, default="points.min")
    parser.add_argument('--Pointsts','-q', help = 'Name of points.ts file', type = str, default="points.ts")
    parser.add_argument('--endianness', help = 'set the endianness of the binary data.  Can be "<" for little-endian or ">" for big-endian', type = str, default="=")
    args = parser.parse_args()

    db = Database()
    converter = OptimDBConverter(args.ndof, db, mindata=args.Mindata, 
                 tsdata=args.Tsdata, pointsmin=args.Pointsmin, pointsts=args.Pointsts,
                 endianness=args.endianness, assert_coords=False)
    converter.pointsmin_data = None
    converter.pointsts_data = None
    converter.ReadMindata()
    converter.ReadTSdata()
    
    system = BaseSystem()
    system.get_ndof = lambda : args.ndof
    run_gui(system, db=db)
    
    print db.number_of_minima()


if __name__ == "__main__":
    main()