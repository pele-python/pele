import argparse
import numpy as np

from pele.gui import run_gui
from pele.systems import BaseSystem

from pele.utils.optim_compatibility import OptimDBConverter

desc="""Start a limited functionality pele gui from a PATHSAMPLE database.
The following files must be present: min.data, ts.data, points.min, and points.ts.  
See the PATHSAMPLE documentation for descriptions of each of those files.
"""  

def main():
    parser = argparse.ArgumentParser(description=desc, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--ndof', help='Number of total degrees of freedom (e.g. 3*number of atoms).  This is simply the length of a coordinates vector.', 
                        type=int, default=1)
#    parser.add_argument('--Database','-d', help = 'Name of database to write into', type = str, default="optimdb.sqlite")
#    parser.add_argument('--Mindata','-m', help = 'Name of min.data file', type = str, default="min.data")
#    parser.add_argument('--Tsdata','-t', help = 'Name of ts.data file', type = str, default="ts.data")
#    parser.add_argument('--Pointsmin','-p', help = 'Name of points.min file', type = str, default="points.min")
#    parser.add_argument('--Pointsts','-q', help = 'Name of points.ts file', type = str, default="points.ts")
#    parser.add_argument('--endianness', help = 'set the endianness of the binary data.  Can be "<" for little-endian or ">" for big-endian', type = str, default="=")
    args = parser.parse_args()

    system = BaseSystem()
    db = system.create_database()
    converter = OptimDBConverter(db, mindata="min.data", 
                 tsdata="ts.data", pointsmin="points.min", pointsts="points.ts",
                 endianness="=", assert_coords=False)
    converter.convert_no_coords()
    
    system.get_ndof = lambda : args.ndof
    run_gui(system, db=db)
    
    print db.number_of_minima()


if __name__ == "__main__":
    main()