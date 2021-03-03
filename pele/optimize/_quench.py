"""
wrappers for the various optimizers.

warning: we've tried to make the most common options the same but there are still differences 

we should make this consistent with scipy.
scipy.minimize would do a similar thing
"""
from __future__ import print_function

import numpy as np

from pele.optimize import LBFGS, MYLBFGS, Fire, Result, LBFGS_CPP, ModifiedFireCPP

__all__ = ["lbfgs_scipy", "fire", "lbfgs_py", "mylbfgs", "cg",
           "steepest_descent", "bfgs_scipy", "lbfgs_cpp"]


def lbfgs_scipy(coords, pot, iprint=-1, tol=1e-3, nsteps=15000):
    """
    a wrapper function for lbfgs routine in scipy
    
    .. warn::
        the scipy version of lbfgs uses linesearch based only on energy
        which can make the minimization stop early.  When the step size
        is so small that the energy doesn't change to within machine precision (times the
        parameter `factr`) the routine declares success and stops.  This sounds fine, but
        if the gradient is analytical the gradient can still be not converged.  This is
        because in the vicinity of the minimum the gradient changes much more rapidly then
        the energy.  Thus we want to make factr as small as possible.  Unfortunately,
        if we make it too small the routine realizes that the linesearch routine
        isn't working and declares failure and exits.
        
        So long story short, if your tolerance is very small (< 1e-6) this routine
        will probably stop before truly reaching that tolerance.  If you reduce `factr` 
        too much to mitigate this lbfgs will stop anyway, but declare failure misleadingly.  
    """
    import scipy.optimize

    res = Result()
    res.coords, res.energy, dictionary = scipy.optimize.fmin_l_bfgs_b(pot.getEnergyGradient,
                                                                      coords, iprint=iprint, pgtol=tol, maxfun=nsteps,
                                                                      factr=10.)
    res.grad = dictionary["grad"]
    res.nfev = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    # res.nsteps = dictionary['nit'] #  new in scipy version 0.12
    res.nsteps = res.nfev
    res.message = dictionary['task']
    res.success = True
    if warnflag > 0:
        print("warning: problem with quench: ", end=' ')
        res.success = False
        if warnflag == 1:
            res.message = "too many function evaluations"
        else:
            res.message = str(dictionary['task'])
        print(res.message)
    # note: if the linesearch fails the lbfgs may fail without setting warnflag.  Check
    # tolerance exactly
    if False:
        if res.success:
            maxV = np.max(np.abs(res.grad))
            if maxV > tol:
                print("warning: gradient seems too large", maxV, "tol =", tol, ". This is a known, but not understood issue of scipy_lbfgs")
                print(res.message)
    res.rms = res.grad.std()
    return res


def fire(coords, pot, tol=1e-3, nsteps=100000, **kwargs):
    """
    A wrapper function for the pele FIRE implementation
    """
    opt = Fire(coords, pot, **kwargs)
    res = opt.run(fmax=tol, steps=nsteps)
    return res


def cg(coords, pot, iprint=-1, tol=1e-3, nsteps=5000, **kwargs):
    """
    a wrapper function for conjugate gradient routine in scipy
    """
    import scipy.optimize

    ret = scipy.optimize.fmin_cg(pot.getEnergy, coords, pot.getGradient,
                                 gtol=tol, full_output=True, disp=iprint > 0,
                                 maxiter=nsteps, **kwargs)
    res = Result()
    res.coords = ret[0]
    res.nfev = ret[2]
    res.nfev += ret[3]  # calls to gradient
    res.success = True
    warnflag = ret[4]
    if warnflag > 0:
        # print "warning: problem with quench: ",
        res.success = False
        if warnflag == 1:
            res.message = "Maximum number of iterations exceeded"
        if warnflag == 2:
            print("Gradient and/or function calls not changing")
    res.energy, res.grad = pot.getEnergyGradient(res.coords)
    res.nfev += 1
    g = res.grad
    res.rms = np.linalg.norm(g) / np.sqrt(len(g))
    return res


def steepest_descent(x0, pot, iprint=-1, dx=1e-4, nsteps=100000,
                     tol=1e-3, maxstep=-1., events=None):
    """steepest descent minimization
    
    Notes
    -----
    this should never be used except for testing purposes.  It is a bad implementation
    of a terrible minimization routine.  It will be very slow.
    """
    N = len(x0)
    x = x0.copy()
    E, V = pot.getEnergyGradient(x)
    funcalls = 1
    for k in range(nsteps):
        stp = -V * dx
        if maxstep > 0:
            stpsize = np.max(np.abs(V))
            if stpsize > maxstep:
                stp *= maxstep / stpsize
        x += stp
        E, V = pot.getEnergyGradient(x)
        funcalls += 1
        rms = np.linalg.norm(V) / np.sqrt(N)
        if iprint > 0:
            if funcalls % iprint == 0:
                print("step %8d energy %20.12g rms gradient %20.12g" % (funcalls, E, rms))
        if events is not None:
            for event in events:
                event(energy=E, coords=x, rms=rms)
        if rms < tol:
            break
    res = Result()
    res.coords = x
    res.energy = E
    res.rms = rms
    res.grad = V
    res.nfev = funcalls
    res.nsteps = k
    res.success = res.rms <= tol
    return res


def bfgs_scipy(coords, pot, iprint=-1, tol=1e-3, nsteps=5000, **kwargs):
    """
    a wrapper function for the scipy BFGS algorithm
    """
    import scipy.optimize

    ret = scipy.optimize.fmin_bfgs(pot.getEnergy, coords, fprime=pot.getGradient,
                                   gtol=tol, full_output=True, disp=iprint > 0,
                                   maxiter=nsteps, **kwargs)
    res = Result()
    res.coords = ret[0]
    res.energy = ret[1]
    res.grad = ret[2]
    res.rms = np.linalg.norm(res.grad) / np.sqrt(len(res.grad))
    res.nfev = ret[4] + ret[5]
    res.nsteps = res.nfev  # not correct, but no better information
    res.success = np.max(np.abs(res.grad)) < tol
    return res


def lbfgs_py(coords, pot, **kwargs):
    lbfgs = LBFGS(coords, pot, **kwargs)
    return lbfgs.run()


def lbfgs_cpp(coords, pot, **kwargs):
    lbfgs = LBFGS_CPP(coords, pot, **kwargs)
    return lbfgs.run()


def mylbfgs(coords, pot, **kwargs):
    lbfgs = MYLBFGS(coords, pot, **kwargs)
    return lbfgs.run()


def modifiedfire_cpp(coords, pot, **kwargs):
    modifiedfire = ModifiedFireCPP(coords, pot, **kwargs)
    return modifiedfire.run()

