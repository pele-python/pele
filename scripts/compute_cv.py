"""
a script to calculate Cv from the Harmonic Superposition Approximation
"""

import argparse
import numpy as np
from pygmin.thermodynamics.heat_capacity import minima_to_cv
from pygmin.storage import Database


        
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
    args = parser.parse_args()
    print args.fname
    k = args.k
    
    dbfname = args.fname
    db = Database(dbfname)

    Tmin = args.Tmin
    Tmax = args.Tmax
    nT = args.Tcount
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    Z, U, U2, Cv = minima_to_cv(db.minima(), T, k)
    
    with open(args.o, "w") as fout:
        fout.write("#T Cv <E> <E**2>\n")
        for vals in zip(T, Cv, U, U2):
            fout.write("%g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.plot(T, Cv, '-')
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig(args.o + ".pdf")
        