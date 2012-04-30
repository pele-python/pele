import scipy.optimize.lbfgsb

def quench(coords, getEnergyGradient, iprint = -1, tol = 1e-3):
    """
    a wrapper function for lbfgs routine in scipy
    """
    newcoords, newE, dictionary = scipy.optimize.fmin_l_bfgs_b(getEnergyGradient, coords, iprint=iprint, pgtol=tol)
    V = dictionary["grad"]
    funcalls = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    if warnflag > 0:
        print "warning: problem with quench: ",
        if warnflag == 1:
            print "too many function evaluations"
        else:
            print dictionary['task']
    rms = V.std()
    return newcoords, newE, rms, funcalls 

def fire(coords, getEnergyGradient, iprint = -1, tol = 1e-3):
    import optimize.fire as fire
    import math
    import numpy as np
    opt = fire.Fire(coords, getEnergyGradient)
    opt.run()
    e,g = getEnergyGradient(opt.coords)
    rms = np.linalg.norm(g)/math.sqrt(len(g))
    return opt.coords, e, rms, opt.nsteps
