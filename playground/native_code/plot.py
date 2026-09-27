from past.utils import old_div
import numpy as np
import pylab as pl

x = np.loadtxt("results.txt").transpose()
pl.plot(x[0], old_div(x[2], x[1]), "o-")
pl.ylabel("time old / time new")
pl.xlabel("number of atoms")
pl.title("100 quenches lj system")
pl.gca().set_yscale("log")
pl.show()
