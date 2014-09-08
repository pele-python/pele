"""
a script to calculate Cv from the Harmonic Superposition Approximation
"""

import argparse
import numpy as np
from pele.thermodynamics import minima_to_cv
from pele.storage import Database
from pele.utils.optim_compatibility import OptimDBConverter

import matplotlib as mpl
mpl.use("Agg") # so we can use it without an X server
import matplotlib.pyplot as plt


        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute the heat capacity from a database of minima using the harmonic" 
                                     "superposition approximation.  The output will be written to Cv and Cv.pdf unless"
                                     "specified differently")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("fname", type=str, help="Database file name")
    parser.add_argument("-k", type=int, help="Number of vibrational degrees freedom.  "
                        "If this is not specified the heat capacity will be off by an additive constant", default=0)
    parser.add_argument("-o", metavar="output", type=str, 
                        help="output file name. The heat capacity will be written to output and output.pdf", default="Cv")   
    parser.add_argument("--Tmin", type=float, help="Minimum temperature for the calculation.", default=0.02)
    parser.add_argument("--Tmax", type=float, help="Minimum temperature for the calculation.", default=1.0)
    parser.add_argument("--Tcount", type=int, help="Number of temperature points for the calculation.", default=300)
    parser.add_argument("--OPTIM", action="store_true", help="read data from a min.data file instead."
                        "fname should be the filename of the min.data file")
    args = parser.parse_args()
    print args.fname
    print args
    k = args.k
    
    # get the list of minima
    if args.OPTIM:
        # fname is a min.data file
        db = Database()
        converter = OptimDBConverter(db, mindata=args.fname)
        converter.convert_no_coords()
        minima = db.minima()
    else:
        dbfname = args.fname
        db = Database(dbfname, createdb=False)
        minima = [m for m in db.minima() if m.fvib is not None and m.pgorder is not None]
        if len(minima) == 0:
            print "There are not minima with the necessary thermodynamic information in the database.  Have you computed the normal mode"\
                  " frequencies and point group order for all the minima?  See pele.thermodynamics "\
                  " for more information"
            exit(1)
    print "computing heat capacity from", len(minima), "minima"

    Tmin = args.Tmin
    Tmax = args.Tmax
    nT = args.Tcount
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    Z, U, U2, Cv = minima_to_cv(minima, T, k)
    
    with open(args.o, "w") as fout:
        fout.write("#T Cv <E> <E**2>\n")
        for vals in zip(T, Cv, U, U2):
            fout.write("%g %g %g %g\n" % vals)
    
    plt.plot(T, Cv, '-')
    plt.xlabel("T")
    plt.ylabel("Cv")
    plt.savefig(args.o + ".pdf")
        
