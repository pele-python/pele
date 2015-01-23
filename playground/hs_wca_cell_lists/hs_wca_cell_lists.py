from __future__ import division
import numpy as np
from pele.potentials import HS_WCA
from pele.optimize import ModifiedFireCPP, LBFGS_CPP
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib 
from matplotlib.patches import Circle 
import pylab 
from scipy.optimize import curve_fit
import copy
from pele.systems import put_in_box

use_periodic = True
use_periodic_frozen = False

def save_pdf(plt, file_name):
    pdf = PdfPages(file_name)
    plt.savefig(pdf, format="pdf")
    pdf.close()
    plt.close()
    
def to_string(inp, digits_after_point = 16):
    format_string = "{0:."
    format_string += str(digits_after_point)
    format_string += "f}"
    return format_string.format(inp)

def plot_disks(coords, radii, boxv, colors=None, sca=0):
    def myscatter(ax, colormap, x, y, radii, colors): 
        for x1,y1,r,c in zip(x, y, radii, colormap(colors)): 
            ax.add_patch(Circle((x1,y1), r, fc=c)) 
    put_in_box(coords, boxv)
    coords = np.reshape(coords, (len(radii),2))
    fig=pylab.figure() 
    ax=fig.add_subplot(111, aspect='equal') 
    myscatter(ax, matplotlib.cm.jet, coords[:,0],coords[:,1], radii*sca, np.ones(len(radii)))
    myscatter(ax, matplotlib.cm.jet, coords[:,0],coords[:,1], radii, np.ones(len(radii))*-1) 
    ax.axis('equal')
    if boxv is not None:
        ax.axes.set_xlim([-boxv[0]/2,boxv[0]/2])
        ax.axes.set_ylim([-boxv[1]/2,boxv[1]/2])
        ax.plot([-boxv[0]/2, -boxv[0]/2], [-boxv[1]/2,boxv[1]/2], color='k', linestyle='-', linewidth=2)
        ax.plot([boxv[0]/2, boxv[0]/2], [-boxv[1]/2,boxv[1]/2], color='k', linestyle='-', linewidth=2)
        ax.plot([-boxv[0]/2, boxv[0]/2], [-boxv[1]/2,-boxv[1]/2], color='k', linestyle='-', linewidth=2)
        ax.plot([-boxv[0]/2, boxv[0]/2], [boxv[1]/2,boxv[1]/2], color='k', linestyle='-', linewidth=2)
    pylab.show()

def reduce_coordinates(mylist, indexes, bdim):
    newlist = mylist.copy().tolist()
    for index in sorted(indexes, reverse=True):
        del newlist[index*bdim:index*bdim+bdim]
    return np.array(newlist)
     
class Config2D(object):
    def __init__(self, nparticles_x, amplitude):
        self.ndim = 2
        self.LX = nparticles_x
        self.LY = self.LX
        self.nparticles_x = nparticles_x
        self.N = self.nparticles_x ** self.ndim
        self.dof = self.ndim * self.N
        self.amplitude = amplitude
        self.x = np.zeros(self.dof)
        for particle in xrange(self.N):
            pid = self.ndim * particle
            self.x[pid] = particle % self.LX
            self.x[pid + 1] = int(particle / self.LX)
        self.x_initial = np.asarray([xi + np.random.uniform(- self.amplitude, self.amplitude) for xi in self.x])
        self.x_initial = np.reshape(self.x_initial, (self.N,2))
        self.x_initial[:,0] -= np.mean(self.x_initial[:,0])
        self.x_initial[:,1] -= np.mean(self.x_initial[:,1])
        self.x_initial = self.x_initial.flatten()  
        #self.radius = 0.3
        #self.sca = 1.5
        self.radius = 0.25
        self.sca = 1.8
        self.radii = np.ones(self.N) * self.radius
        self.eps = 1.0
        self.boxvec = np.array([self.LX, self.LY])
        self.potential = HS_WCA(use_periodic=use_periodic, eps=self.eps,
                         sca=self.sca, radii=self.radii.copy(), ndim=self.ndim, boxvec=self.boxvec.copy())
        self.potential_ = HS_WCA(use_periodic=use_periodic, eps=self.eps,
                         sca=self.sca, radii=self.radii.copy(), ndim=self.ndim, boxvec=self.boxvec.copy())
        self.rcut = 2 * (1 + self.sca) * self.radius
        self.ncellx_scale = 1
        self.potential_cells = HS_WCA(use_periodic=use_periodic,
                               use_cell_lists=True, eps=self.eps,
                               sca=self.sca, radii=self.radii.copy(),
                               boxvec=self.boxvec.copy(),
                               rcut=self.rcut, ndim=self.ndim,
                               ncellx_scale=self.ncellx_scale)
        self.potential_cells_ = HS_WCA(use_periodic=use_periodic,
                                use_cell_lists=True, eps=self.eps,
                                sca=self.sca, radii=self.radii.copy(),
                                boxvec=self.boxvec.copy(),
                                rcut=self.rcut, ndim=self.ndim,
                                ncellx_scale=self.ncellx_scale)
        self.tol = 1e-7
        self.maxstep = np.amax(self.radii)
        self.nstepsmax = int(1e6)
        assert(self.boxvec[0]==self.boxvec[1])
        print "x_initial energy:", self.potential.getEnergy(self.x_initial)
        print "x_initial cells energy:", self.potential_cells.getEnergy(self.x_initial)
        assert(self.potential.getEnergy(self.x_initial) == self.potential_.getEnergy(self.x_initial))
        assert(self.potential_cells.getEnergy(self.x_initial) == self.potential_cells_.getEnergy(self.x_initial))
        #assert abs(self.potential.getEnergy(self.x_initial) - self.potential_cells.getEnergy(self.x_initial)) < 1e-10
        assert np.allclose(self.potential.getEnergy(self.x_initial), self.potential_cells.getEnergy(self.x_initial), rtol=1e-10)
        print self.boxvec
        #plot_disks(self.x_initial, self.radii, self.boxvec, sca=self.sca)
        
    def optimize(self, nr_samples = 1):
        self.optimizer =  ModifiedFireCPP(self.x_initial.copy(), self.potential,
                                         dtmax=1, maxstep=self.maxstep,
                                         tol=self.tol, nsteps=1e8, verbosity=-1, iprint=-1)
        self.optimizer_ = LBFGS_CPP(self.x_initial.copy(), self.potential_)
        self.optimizer_cells = ModifiedFireCPP(self.x_initial.copy(), self.potential_cells,
                                         dtmax=1, maxstep=self.maxstep,
                                         tol=self.tol, nsteps=1e8, verbosity=-1, iprint=-1)
        self.optimizer_cells_ = LBFGS_CPP(self.x_initial.copy(), self.potential_cells_)
        print "initial E, x:", self.potential.getEnergy(self.x_initial.copy())
        print "initial E, x_:", self.potential_cells.getEnergy(self.x_initial.copy())
        
        t0 = time.time()
        print "self.optimizer.run(self.nstepsmax)", self.nstepsmax
        self.optimizer.run(self.nstepsmax)
        self.res_x_final = self.optimizer.get_result()
        t1 = time.time()
        self.optimizer_cells.run(self.nstepsmax)
        self.res_x_final_cells = self.optimizer_cells.get_result()
        t2 = time.time()
                
        self.x_final = self.res_x_final.coords
        self.x_final_cells = self.res_x_final_cells.coords
        print "fire final E, x:", self.optimizer.get_result().energy
        print "fire final E, x_cells:", self.optimizer_cells.get_result().energy
        print "fire final E, plain: ", self.potential.getEnergy(self.x_final)
        print "fire final E, cell: ", self.potential_cells.getEnergy(self.x_final_cells)
        print "fire number of particles:", self.N
        print "fire time no cell lists:", t1 - t0, "sec"
        print "fire time cell lists:", t2 - t1, "sec"
        print "fire ratio:", (t1 - t0) / (t2 - t1)
        
        if not self.res_x_final.success or not self.res_x_final_cells.success:
            print "-------------"
            print "res_x_final.rms:", self.res_x_final.rms
            print "res_x_final.nfev:", self.res_x_final.nfev
            print "res_x_final_cells.rms:", self.res_x_final_cells.rms
            print "res_x_final_cells.nfev:", self.res_x_final_cells.nfev
            print "self.res_x_final.success", self.res_x_final.success
            print "self.res_x_final_cells.success", self.res_x_final_cells.success
            print "-------------"
            plot_disks(self.x_initial, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_cells, self.radii, self.boxvec, sca=self.sca)
                
        self.optimizer_.run(self.nstepsmax)
        self.res_x_final_ = self.optimizer_.get_result()
        t3 = time.time()
        self.optimizer_cells_.run(self.nstepsmax)
        self.res_x_final_cells_ = self.optimizer_cells_.get_result()
        t4 = time.time()
        
        self.x_final_ = self.res_x_final_.coords
        self.x_final_cells_ = self.res_x_final_cells_.coords
        print "lbfgs final E, x:", self.optimizer_.get_result().energy
        print "lbfgs final E, x_cells:", self.optimizer_cells_.get_result().energy
        print "lbfgs final E, plain: ", self.potential_.getEnergy(self.x_final_)
        print "lbfgs final E, cell: ", self.potential_cells_.getEnergy(self.x_final_cells_)
        print "lbfgs number of particles:", self.N
        print "lbfgs time no cell lists:", t3 - t2, "sec"
        print "lbfgs time cell lists:", t4 - t3, "sec"
        print "lbfgs ratio:", (t3 - t2) / (t4 - t3)
        
        if not self.res_x_final_.success or not self.res_x_final_cells_.success or not self.res_x_final.success or not self.res_x_final_cells.success:
            print "-------------"
            print "res_x_final_.rms:", self.res_x_final_.rms
            print "res_x_final_.nfev:", self.res_x_final_.nfev
            print "res_x_final_cells_.rms:", self.res_x_final_cells_.rms
            print "res_x_final_cells_.nfev:", self.res_x_final_cells_.nfev
            print "self.res_x_final_.success", self.res_x_final.success
            print "self.res_x_final_cells_.success", self.res_x_final_cells.success
            print "-------------"
            plot_disks(self.x_initial, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_cells_, self.radii, self.boxvec, sca=self.sca)
        
        assert(self.res_x_final.success)
        assert(self.res_x_final_cells.success)
        assert(self.res_x_final_.success)
        assert(self.res_x_final_cells_.success)
        
        for (xci, xi) in zip(self.x_final_cells, self.x_final):
            passed = (np.abs(xci - xi) < 1e-10)
            if (passed is False):
                print "xci", xci
                print "xi", xi
                assert(passed)
        print "energy no cell lists:", self.res_x_final.energy
        print "energy cell lists:", self.res_x_final_cells.energy
        self.t_ratio = (t1 - t0) / (t2 - t1)
        self.t_ratio_lbfgs = (t3 - t2) / (t4 - t3)
        
class Config2DFrozenBoundary(object):
    def __init__(self, nparticles_x, amplitude):
        self.ndim = 2
        self.LX = nparticles_x
        self.LY = self.LX
        self.nparticles_x = nparticles_x
        self.N = self.nparticles_x ** self.ndim
        self.amplitude = amplitude
        self.dof = self.ndim * self.N
        self.x = np.zeros(self.dof)
        self.frozen_atoms = []
        for particle in xrange(self.N):
            pid = self.ndim * particle
            xcoor = particle % self.LX
            ycoor = int(particle / self.LX)
            self.x[pid] = xcoor
            self.x[pid + 1] = ycoor
            if xcoor == 0 or xcoor == self.LX - 1 or ycoor == 0 or ycoor == self.LY - 1:
                self.frozen_atoms.append(particle)
        self.x_initial = copy.copy(self.x)
        for particle in xrange(self.N):
            if particle not in self.frozen_atoms:
                pid = self.ndim * particle
                self.x_initial[pid] += np.random.uniform(- self.amplitude, self.amplitude)
                self.x_initial[pid + 1] += np.random.uniform(- self.amplitude, self.amplitude)
        self.x_initial = np.reshape(self.x_initial, (self.N,2))
        self.x_initial[:,0] -= np.mean(self.x_initial[:,0])
        self.x_initial[:,1] -= np.mean(self.x_initial[:,1])
        self.x_initial = self.x_initial.flatten()
        min_x = np.amin(self.x_initial)
        if min_x < 0:
            self.x_initial -= min_x
        #self.radius = 0.3
        #self.sca = 1.5
        self.radius = 0.25
        self.sca = 1.8
        self.radii = np.ones(self.N) * self.radius
        self.eps = 1.0
        max_edge = np.amax([np.amax(self.x_initial), np.abs(np.amin(self.x_initial))]) + 2 * self.amplitude + (1 + self.sca) * self.radius
        self.boxvec = np.array([max_edge, max_edge])
        self.frozen_atoms1 = np.array(self.frozen_atoms)
        self.frozen_atoms2 = np.array(self.frozen_atoms)
        print "self.frozen_atoms1", self.frozen_atoms1
        self.potential = HS_WCA(use_frozen=True, use_periodic=use_periodic_frozen,
                         reference_coords=self.x_initial,
                         frozen_atoms=self.frozen_atoms1,
                         eps=self.eps, sca=self.sca, radii=self.radii,
                         ndim=self.ndim, boxvec=self.boxvec)
        self.rcut =  2 * (1 + self.sca) * self.radius
        self.ncellx_scale = 1.0
        self.potential_cells = HS_WCA(use_frozen=True,
                               use_periodic=use_periodic_frozen, use_cell_lists=True,
                               eps=self.eps, sca=self.sca,
                               radii=self.radii, boxvec=self.boxvec,
                               reference_coords=self.x_initial,
                               rcut=self.rcut, ndim=self.ndim,
                               ncellx_scale=self.ncellx_scale,
                               frozen_atoms=self.frozen_atoms2)
        self.tol = 1e-7
        self.maxstep = np.amax(self.radii)
        self.nstepsmax = int(1e6)
        assert(self.boxvec[0]==self.boxvec[1])
        self.x_initial_red = reduce_coordinates(self.x_initial,self.frozen_atoms,self.ndim)
        print "x_initial energy:", self.potential.getEnergy(self.x_initial_red)
        print "x_initial cells energy:", self.potential_cells.getEnergy(self.x_initial_red)
        #assert abs(self.potential.getEnergy(self.x_initial_red) - self.potential_cells.getEnergy(self.x_initial_red)) < 1e-10
        assert np.allclose(self.potential.getEnergy(self.x_initial_red), self.potential_cells.getEnergy(self.x_initial_red), rtol=1e-10)
        print self.boxvec
    
    def optimize(self, nr_samples=1):
        self.x_initial_red = reduce_coordinates(self.x_initial,self.frozen_atoms,self.ndim)
        self.optimizer = ModifiedFireCPP(self.x_initial_red.copy(), self.potential, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_ = LBFGS_CPP(self.x_initial_red.copy(), self.potential)
        self.optimizer_cells = ModifiedFireCPP(self.x_initial_red.copy(), self.potential_cells, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_cells_ = LBFGS_CPP(self.x_initial_red.copy(), self.potential_cells)
        t0 = time.time()
        print "self.optimizer.run(self.nstepsmax)", self.nstepsmax
        self.optimizer.run(self.nstepsmax)
        self.res_x_final = self.optimizer.get_result()
        t1 = time.time()
        self.optimizer_cells.run(self.nstepsmax)
        self.res_x_final_cells = self.optimizer_cells.get_result()
        t2 = time.time()
                
        self.x_final = self.res_x_final.coords
        self.x_final_cells = self.res_x_final_cells.coords
        print "fire final E, x:", self.optimizer.get_result().energy
        print "fire final E, x_cells:", self.optimizer_cells.get_result().energy
        print "fire final E, plain: ", self.potential.getEnergy(self.x_final)
        print "fire final E, cell: ", self.potential_cells.getEnergy(self.x_final_cells)
        print "fire number of particles:", self.N
        print "fire time no cell lists:", t1 - t0, "sec"
        print "fire time cell lists:", t2 - t1, "sec"
        print "fire ratio:", (t1 - t0) / (t2 - t1)
        
        if not self.res_x_final.success or not self.res_x_final_cells.success:
            print "-------------"
            print "res_x_final.rms:", self.res_x_final.rms
            print "res_x_final.nfev:", self.res_x_final.nfev
            print "res_x_final_cells.rms:", self.res_x_final_cells.rms
            print "res_x_final_cells.nfev:", self.res_x_final_cells.nfev
            print "self.res_x_final.success", self.res_x_final.success
            print "self.res_x_final_cells.success", self.res_x_final_cells.success
            print "-------------"
            plot_disks(self.x_initial, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_cells, self.radii, self.boxvec, sca=self.sca)
                
        self.optimizer_.run(self.nstepsmax)
        self.res_x_final_ = self.optimizer_.get_result()
        t3 = time.time()
        self.optimizer_cells_.run(self.nstepsmax)
        self.res_x_final_cells_ = self.optimizer_cells_.get_result()
        t4 = time.time()
        
        self.x_final_ = self.res_x_final_.coords
        self.x_final_cells_ = self.res_x_final_cells_.coords
        print "lbfgs final E, x:", self.optimizer_.get_result().energy
        print "lbfgs final E, x_cells:", self.optimizer_cells_.get_result().energy
        print "lbfgs final E, plain: ", self.potential.getEnergy(self.x_final_)
        print "lbfgs final E, cell: ", self.potential_cells.getEnergy(self.x_final_cells_)
        print "lbfgs number of particles:", self.N
        print "lbfgs time no cell lists:", t3 - t2, "sec"
        print "lbfgs time cell lists:", t4 - t3, "sec"
        print "lbfgs ratio:", (t3 - t2) / (t4 - t3)
        
        if not self.res_x_final_.success or not self.res_x_final_cells_.success or not self.res_x_final.success or not self.res_x_final_cells.success:
            print "-------------"
            print "res_x_final_.rms:", self.res_x_final_.rms
            print "res_x_final_.nfev:", self.res_x_final_.nfev
            print "res_x_final_cells_.rms:", self.res_x_final_cells_.rms
            print "res_x_final_cells_.nfev:", self.res_x_final_cells_.nfev
            print "self.res_x_final_.success", self.res_x_final.success
            print "self.res_x_final_cells_.success", self.res_x_final_cells.success
            print "-------------"
            plot_disks(self.x_initial, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_, self.radii, self.boxvec, sca=self.sca)
            plot_disks(self.x_final_cells_, self.radii, self.boxvec, sca=self.sca)
        
        assert(self.res_x_final.success)
        assert(self.res_x_final_cells.success)
        assert(self.res_x_final_.success)
        assert(self.res_x_final_cells_.success)
        for (xci, xi) in zip(self.x_final_cells, self.x_final):
            passed = (np.abs(xci - xi) < 1e-10)
            if (passed is False):
                print "xci", xci
                print "xi", xi
                assert(passed)
        print "energy no cell lists:", self.res_x_final.energy
        print "energy cell lists:", self.res_x_final_cells.energy
        self.t_ratio = (t1 - t0) / (t2 - t1)
        self.t_ratio_lbfgs = (t3 - t2) / (t4 - t3)

def ll(x, a, b):
    return a * x ** b

def measurement(nr_samples, LXmax, amplitude):
    nr_particles = []
    time_ratio = []
    time_ratio_lbfgs = []
    for LX in xrange(3, LXmax + 1):
        print "---LX: ", LX
        t_fire = 0
        t_lbfgs = 0
        nrp = 0
        for tmp in xrange(nr_samples):
            conf = Config2D(LX, amplitude)
            conf.optimize()
            t_fire = t_fire + conf.t_ratio / nr_samples
            t_lbfgs = t_lbfgs + conf.t_ratio_lbfgs / nr_samples
            nrp = conf.N
        nr_particles.append(nrp)
        time_ratio.append(t_fire)
        time_ratio_lbfgs.append(t_lbfgs)    
        print "---done:", LX, "max: ", LXmax  
    popt, pcov = curve_fit(ll, nr_particles, time_ratio)
    plt.plot(nr_particles, np.array(time_ratio) - 1, "s", label = r"M-FIRE")
    plt.plot(nr_particles, np.array(time_ratio_lbfgs) - 1, "o", label = r"LBFGS")
    xr = np.linspace(nr_particles[0], nr_particles[-1])
    plt.plot(xr, np.array([ll(xi, popt[0], popt[1]) for xi in xr]) - 1, "k", label=r"Exponent: " + to_string(popt[1], 2))
    plt.axhline(0, color='black')
    ax = plt.axes()
    ax.grid(True)
    plt.legend(loc = 2)
    plt.xlabel(r"Number of particles")
    plt.ylabel(r"Time no cell lists / time cell lists - 1")
    save_pdf(plt, "hs_wca_cell_lists.pdf")
    print "done: non-frozen measurement"

def measurement_frozen(nr_samples, LXmax, amplitude):
    nr_particles = []
    time_ratio = []
    time_ratio_lbfgs = []
    for LX in xrange(3, LXmax + 1):
        print "---LX: ", LX
        t_fire = 0
        t_lbfgs = 0
        nrp = 0
        for tmp in xrange(nr_samples):
            print "frozen sample", tmp
            conf = Config2DFrozenBoundary(LX, amplitude)
            print "construction"
            conf.optimize()
            print "optimization"
            t_fire = t_fire + conf.t_ratio / nr_samples
            t_lbfgs = t_lbfgs + conf.t_ratio_lbfgs / nr_samples
            nrp = conf.N
        nr_particles.append(nrp)
        time_ratio.append(t_fire)
        time_ratio_lbfgs.append(t_lbfgs)    
        print "---done:", LX, "max: ", LXmax  
    popt, pcov = curve_fit(ll, nr_particles, time_ratio)
    plt.plot(nr_particles, np.array(time_ratio)-1, "s", label = r"M-FIRE")
    plt.plot(nr_particles, np.array(time_ratio_lbfgs)-1, "o", label = r"LBFGS")
    xr = np.linspace(nr_particles[0], nr_particles[-1])
    plt.plot(xr, np.array([ll(xi, popt[0], popt[1]) for xi in xr])-1, "k", label=r"Exponent: " + to_string(popt[1], 2))
    plt.axhline(0, color='black')
    ax = plt.axes()
    ax.grid(True)
    plt.legend(loc = 2)
    plt.xlabel(r"Number of particles")
    plt.ylabel(r"Time no cell lists / time cell lists - 1")
    save_pdf(plt, "hs_wca_cell_lists_frozen.pdf")
    print "done: frozen measurement"

if __name__ == "__main__":
    
    np.random.seed(42)
    measurement(20, 14, 0.1)
    measurement_frozen(20, 14, 0.1)
