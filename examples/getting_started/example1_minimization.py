"""

First steps
-----------
`pele` is first and foremost a library.  There is no single executable that you
can call.  Instead, you use python as a scripting language to explicitly state
what you want to do, importing the tools from `pele` to do it.

As a very simple example, let's generate a random Lennard-Jones configuration and
minimize it.
::
    import numpy as np
    from pele.potentials import LJ
    from pele.optimize import lbfgs_py
    
    natoms = 5
    x = np.random.uniform(-2, 2, natoms*3)
    pot = LJ()
    result = lbfgs_py(x, pot)
    print result

In the above we used three imports.  We used the `numpy` library to construct a random
one dimensional array.  We used the Lennard-Jones potential :class:`.LJ`, and we used the minimization
routine :func:`.lbfgs_py` which is just a wrapper for the class :class:`.LBFGS`.

The return value is an optimization result, which is just a container (:class:`.Result`) that stores
the final energy, final coordinates, the number of function calls, etc.

If we want to then save the minimized coordinates in an xyz file we can use the function :func:`.write_xyz`
::
    from pele.utils.xyz import write_xyz
    with open("out.xyz", "w") as fout:
        title = "energy = " + str(result.energy)
        write_xyz(fout, result.coords, title=title)

"""
import numpy as np
from pele.potentials import LJ
from pele.optimize import lbfgs_py

natoms = 5
x = np.random.uniform(-2, 2, natoms*3)
pot = LJ()
result = lbfgs_py(x, pot)
print result


from pele.utils.xyz import write_xyz
with open("out.xyz", "w") as fout:
    title = "energy = " + str(result.energy)
    write_xyz(fout, result.coords, title=title)
