import numpy as np
from pygmin.optimize import quench

xt=[]
tt=[]

# add event callback
def findTS(x0, potential, direction=None, tol=1.0e-6, maxstep=0.1, **kwargs):
    '''
    Wrapper for DimerSearch to find transition states
    '''    
    search = DimerSearch(x0, potential, direction=direction, **kwargs)
    x, E, rms, tmp = quench.mylbfgs(x0, search.getEnergyGradient, tol=tol, maxstep=maxstep, maxErise=10.)
    
    ret = dict()
    ret["eigenvec"]=search.tau
    # TODO
    ret["eigenval"]=0.0
    ret["rms"]=rms
    
    return x,E,ret

class DimerSearch(object):
    '''
    single ended transition search method using hybrid eigenvector following / dimer method. 
    This implementation follows 
    "Superlinearly converging dimer method for transition state search"
    Johannes Kaestner and Paul Sherwood, J. Chem. Phys. 128, 014106 (2008)
    http://dx.doi.org/10.1063/1.2815812
    '''


    def __init__(self, x0, potential, direction=None, delta=1e-5, max_rotsteps=100, theta_cut=0.05, zeroEigenVecs=None):
        '''
        Constructor
        '''
        # potential to work with
        self.potential = potential
        # point to start search from
        self.x0 = x0
        # size of the dimer
        self.delta = 1e-3
        # maximum number of dimer rotation updates per step
        self.max_rotsteps = max_rotsteps
        # threshold to perform dimer rotation updates
        self.theta_cut = theta_cut
        # callback to calculate zero eigenvectors
        self.zeroEigenVecs = zeroEigenVecs
        # searches starting in this direction have already been performed
        self.tau_done=[]
        # current list of eigenvectors to projected out
        self.tau_ignore=[]
        
        if(direction==None):
            direction=np.random.random(x0.shape) - 0.5
        self.findNextTS(direction)
    
    def findNextTS(self, direction):
        self.tau = direction/np.linalg.norm(direction)
        E,g = self.potential.getEnergyGradient(self.x0)
        self.updateRotation(self.x0, E, g)
        import copy
        self.tau_ignore=copy.copy(self.tau_done)
        self.tau_done.append([self.tau.copy(), 0.])
        
    def getEnergyGradient(self, x0):
        E,g = self.potential.getEnergyGradient(x0)
        self.updateRotation(x0, E, g)
        
        g = g - 2.*np.dot(g, self.tau)*self.tau
        xt.append(x0)
        tt.append(self.tau)
        return E,g
    
    # should this be iterative?, go to utility with ljsyszem zero eigenvecs
    def orthogonalize(self, x, vecs, vecs2=None):
        for v in vecs:
            x -= np.dot(x, v)*v
        if(vecs2):
            for v in vecs2:
                x -= np.dot(x, v)*v
        
    
    def getOrthogonalGradient(self, x, eigenvecs1, eigenvecs2=None):
        E, g = self.potential.getEnergyGradient(x)
        #self.orthogonalize(g, eigenvecs1, eigenvecs2)
        return g
    
    def updateRotation(self, x0, E0, grad0_):
        iter_rot = 0
        
        # get the zero eigenvalues
        zev = []
        if(self.zeroEigenVecs):
            zev = self.zeroEigenVecs(x0)

        # remove zero eigenvalues from gradient
        grad0 = grad0_.copy()
        self.orthogonalize(grad0, zev, [v[0] for v in self.tau_ignore])

        # update ignore list for eigenvalues
        for t in self.tau_ignore:
            grad1 = self.getOrthogonalGradient(x0 + t[1]*self.delta, zev)
            t[1] = np.dot((grad1 - grad0), t[1])/self.delta
            
        while iter_rot < self.max_rotsteps:
            #self.tau = self.tau/np.linalg.norm(self.tau)
            # construct dimer image and get energy + gradient
            
            # remove zero eigenvalues from tau
            self.orthogonalize(self.tau, zev, [v[0] for v in self.tau_ignore])
            self.tau /= np.linalg.norm(self.tau)
            
            x1 = x0 + self.tau*self.delta
            grad1 = self.getOrthogonalGradient(x1, zev, [v[0] for v in self.tau_ignore])
            
            # calculate the rotational force of dimer
            F_rot = -2.*(grad1 - grad0) + 2.*np.dot(grad1 - grad0, self.tau)*self.tau
            
            # For now just use steepest descent search direction for rotation.
            # Replace this by LBFGS
            Theta = F_rot / np.linalg.norm(F_rot)
            
            # calculate curvature C and derivative of curvature
            C = np.dot((grad1 - grad0), self.tau)/self.delta
            dC = 2.*np.dot((grad1 - grad0), Theta)/self.delta
            #print C,self.tau
            # calculate estimated rotation angle
            theta1=-0.5*np.arctan(dC/(2.*np.abs(C)))
            
            # do we need to rotate or already smaller than cutoff?
            if np.abs(theta1) < self.theta_cut:
                return
            
            # create rotated trial dimer
            taup = self.rotate(self.tau, Theta, theta1)
            # remove zero eigenvalues from tau
            self.orthogonalize(taup, zev, [v[0] for v in self.tau_ignore])
            taup /= np.linalg.norm(taup)
            
            x1p = x0 + self.delta * taup
            
            # get the new energy and gradient at trial conviguration
            grad1p = self.getOrthogonalGradient(x1p, zev, [v[0] for v in self.tau_ignore])

            #grad1p = -grad1p
            # get curvature for trial point
            Cp = np.dot((grad1p - grad0), taup)/self.delta
            
            # calculate optimum rotation angle theta_min and taumin
            b1 = 0.5*dC
            a1 = (C - Cp + b1*np.sin(2.*theta1))/(1. - np.cos(2.*theta1))
            theta_min=0.5*np.arctan(b1/a1)
            self.tau = self.rotate(self.tau, Theta, theta_min)
            
            # remove zero eigenvalues from tau
            self.orthogonalize(self.tau, zev, [v[0] for v in self.tau_ignore])
            self.tau /= np.linalg.norm(self.tau)
            

            if np.abs(theta_min) < self.theta_cut:
                return
            iter_rot+=1
        
    
    # rotate the current direction vector tau
    def rotate(self, tau, Theta, angle):
        #print "before", np.linalg.norm(tau)
        #print "after", np.linalg.norm(tau*np.cos(angle) + Theta*np.sin(angle))
        return tau*np.cos(angle) + Theta*np.sin(angle)

import nebtesting as test

class potwrap(object):
    def __init__(self, potential):
        self.potential = potential
        self.trajectory=[]
    
    def getEnergyGradient(self, x):
        self.trajectory.append(x)
        return self.potential.getEnergyGradient(x)
    
    def getEnergy(self, x):
        return self.potential.getEnergy(x)
    
if __name__ == "__main__":
    import pylab as pl
    x = np.arange(.5, 5., .05)
    y = np.arange(.5, 5., .05)
    z = np.zeros([len(x), len(y)])
    potential = potwrap(test.leps())
    for i in range(0, len(x)):
        for j in range(0, len(y)):
                z[j, i] = potential.getEnergy([x[i], y[j]])
    x0 = np.array([.75, 1.5]) #np.random.random(3)
    #final = np.array([2., .75]) #np.random.random(3)

    pl.pcolor(x, y, z, vmax=-0.5, cmap=pl.cm.PuBu)
    tau=np.array([0.5,-1.])#np.random.random(2)-0.5
    x1 = x0 + 0.1*tau
    pl.plot([x0[0], x1[0]], [x0[1], x1[1]])
    x0,E,ret = findTS(x0, potential, tau)
    tau=ret["eigenvec"]
    x1 = x0 + 0.1*tau
    pl.plot([x0[0], x1[0]], [x0[1], x1[1]])
    pl.colorbar()
    pl.xlabel("x")
    pl.ylabel("y")
    for i in range(len(xt)):
        pl.plot(xt[i][0], xt[i][1], 'o')
        pl.plot([xt[i][0],xt[i][0] + 0.05*tt[i][0]], [xt[i][1],xt[i][1]+ 0.05*tt[i][1]], '-')
    print E
    import tstools
    potential.trajectory = []
    m1,m2 = tstools.minima_from_ts(potential.getEnergyGradient, x0, tau, displace=1e-2)
    for xt in potential.trajectory:
      pl.plot(xt[0], xt[1], 'x')
    pl.axis(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5)
    
    pl.show()          

        