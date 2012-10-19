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
        parameters for setting up NEB
            important parameters:
                number of images
        parameters for NEB quench
            Note: the potential for the NEB optimization is not Hamiltonian.
                This means that there is no meaningful energy associated with the potential.  
                Therefore, during the optimization, we can do gradient following, 
                but we can't rely on the energy for, e.g. determining step size, 
                which is the default behavior for many optimizers.  This can be worked around by choosing a 
                small step size and a large maxErise, or by using an optimizer that uses only gradients.
                
                scipy.lbfgs_b seems to work with NEB pretty well, but lbfgs_py and mylbfgs tend to fail.  If 
                you must use one of those try, e.g. maxErise = 1., maxstep=0.01, tol=1e-2 
            important parameters:
                nsteps : should be small (roughly 100).  It often will minimize forever without reaching 
                         the required tolerance.  Also, we just want an estimate for the transition
                         state, we must further refine it with transition state search, which is often
                         much faster than NEB.

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

NEBquenchRoutine = quench.quench
NEBquenchParams = dict()
NEBquenchParams["nsteps"] = 100