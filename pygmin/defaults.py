from pygmin.optimize import quench
from pygmin.NEB import dimer

"""
These are an easy way to change the default parameters of the various features of PyGMIN. 

to change the default quenchRoutine to FIRE, simply run
>>> import pygmin.utils.defaults as defaults
>>> from pygmin.optimize.quench import fire
>>> defaults.quenchRoutine = fire

to further change the tolerance for all quenches do

>>> defaults.quenchRoutine["tol"] = 0.01
"""

"""
notes
-----
The following are sections of the code which should have default settings in this way

    normal quenching :
        which routine
        parameters
        
    lowest eigenvector search :
        which routine
            eventually this can become untied to the transition state search
        parameters
            important parameters: 
                specify zero eigenmodes
    
    minimization in the space tangent to lowest eigenvector
        which routine
        parameters
    
    Transition state search:
        which routine
        parameters
            transition state search uses lowest egenvector search and tangent space minimization, so 
            these parameters would be in addition to those

    NEB:
        which routine:
            NEB, DNEB, etc
        parameters
            important parameters:
                nsteps : should be small (100).  It often will minimize forever without reaching 
                         the required tolerance

    double ended connect : 
        parameters:
            This uses NEB, Transition state search, and mindist, so these 
            parameters will be in addition to those.
            important parameters:
          
    single ended connect : 
        parameters

    best alignment between two structures (mindist):
        parameters
            here the parameters would be, e.g. rotational symmetry, translational symmetry,
            permutational symmetry, etc.  In principle the routine is determined by the parameters

"""

quenchRoutine = quench.lbfgs_py
quenchParams = dict()

tsSearchRoutine = dimer
tsSearchParams = dict()

lowestEigenvectorQuenchParams = dict()

NEBquenchRoutine = quench.lbfgs_py
NEBquenchParams = dict()
NEBquenchParams["nsteps"] = 100