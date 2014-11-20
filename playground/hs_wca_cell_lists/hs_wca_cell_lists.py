from __future__ import division
import numpy as np
from pele.potentials import HS_WCA
from pele.optimize import ModifiedFireCPP, LBFGS_CPP
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import copy

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
        #self.radius = 0.3
        #self.sca = 1.5
        self.radius = 0.25
        self.sca = 1.1
        self.radii = np.ones(self.N) * self.radius
        self.eps = 1.0
        self.boxvec = np.array([self.LX, self.LY])
        self.potential = HS_WCA(use_periodic=True, eps=self.eps,
                         sca=self.sca, radii=self.radii, ndim=self.ndim, boxvec=self.boxvec)
        self.potential_ = HS_WCA(use_periodic=True, eps=self.eps,
                         sca=self.sca, radii=self.radii, ndim=self.ndim, boxvec=self.boxvec)
        self.rcut = 2 * (1 + self.sca) * self.radius
        self.ncellx_scale = 1.0
        self.potential_cells = HS_WCA(use_periodic=True,
                               use_cell_lists=True, eps=self.eps,
                               sca=self.sca, radii=self.radii,
                               boxvec=self.boxvec,
                               reference_coords=self.x_initial,
                               rcut=self.rcut, ndim=self.ndim,
                               ncellx_scale=self.ncellx_scale)
        self.potential_cells_ = HS_WCA(use_periodic=True,
                                use_cell_lists=True, eps=self.eps,
                                sca=self.sca, radii=self.radii,
                                boxvec=self.boxvec,
                                reference_coords=self.x_initial,
                                rcut=self.rcut, ndim=self.ndim,
                                ncellx_scale=self.ncellx_scale)
        self.tol = 1e-7
        self.maxstep = 1.0
        self.nstepsmax = int(1e4)
    def optimize(self, nr_samples = 1):
        self.optimizer = ModifiedFireCPP(self.x_initial, self.potential,
                         tol=self.tol, maxstep=self.maxstep)
        self.optimizer_ = LBFGS_CPP(self.x_initial, self.potential_)
        self.optimizer_cells = ModifiedFireCPP(self.x_initial, self.potential_cells,
                         tol=self.tol, maxstep=self.maxstep)
        self.optimizer_cells_ = LBFGS_CPP(self.x_initial, self.potential_cells_)
        print "initial E, x:", self.optimizer.get_result().energy
        print "initial E, x_:", self.optimizer_cells.get_result().energy
        t0 = time.time()
        print "self.optimizer.run(self.nstepsmax)", self.nstepsmax
        self.optimizer.run(self.nstepsmax)
        t1 = time.time()
        self.optimizer_cells.run(self.nstepsmax)
        t2 = time.time()
        self.optimizer_.run(self.nstepsmax)
        t3 = time.time()
        self.optimizer_cells_.run(self.nstepsmax)
        t4 = time.time()
        self.res_x_final = self.optimizer.get_result()
        self.res_x_final_cells = self.optimizer_cells.get_result()
        self.x_final = self.res_x_final.coords
        self.x_final_cells = self.res_x_final_cells.coords
        print "final E, x:", self.optimizer.get_result().energy
        print "final E, x_:", self.optimizer_cells.get_result().energy
        print "number of particles:", self.N
        print "time no cell lists:", t1 - t0, "sec"
        print "time cell lists:", t2 - t1, "sec"
        print "ratio:", (t1 - t0) / (t2 - t1)
        self.res_x_final_ = self.optimizer_.get_result()
        self.res_x_final_cells_ = self.optimizer_cells_.get_result()
        if not self.res_x_final.success or not self.res_x_final_cells.success:
            print "res_x_final.rms:", self.res_x_final.rms
            print "res_x_final.nfev:", self.res_x_final.nfev
            print "res_x_final_cells.rms:", self.res_x_final_cells.rms
            print "res_x_final_cells.nfev:", self.res_x_final_cells.nfev
            print "-------------"
            print "res_x_final_.rms:", self.res_x_final_.rms
            print "res_x_final_.nfev:", self.res_x_final_.nfev
            print "res_x_final_cells_.rms:", self.res_x_final_cells_.rms
            print "res_x_final_cells_.nfev:", self.res_x_final_cells_.nfev
            print "-------------"
            print "self.res_x_final.success", self.res_x_final.success
            print "self.res_x_final_cells.success", self.res_x_final_cells.success
        assert(self.res_x_final.success)
        assert(self.res_x_final_cells.success)
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
        #self.radius = 0.3
        #self.sca = 1.5
        self.radius = 0.25
        self.sca = 1.1
        self.radii = np.ones(self.N) * self.radius
        self.eps = 1.0
        self.boxvec = np.array([self.LX, self.LY])
        self.frozen_atoms1 = np.array(self.frozen_atoms)
        self.frozen_atoms2 = np.array(self.frozen_atoms)
        print "self.frozen_atoms1", self.frozen_atoms1
        self.potential = HS_WCA(use_frozen=True, use_periodic=True,
                         reference_coords=self.x_initial,
                         frozen_atoms=self.frozen_atoms1,
                         eps=self.eps, sca=self.sca, radii=self.radii,
                         ndim=self.ndim, boxvec=self.boxvec)
        self.rcut =  2 * (1 + self.sca) * self.radius
        self.ncellx_scale = 1.0
        self.potential_cells = HS_WCA(use_frozen=True,
                               use_periodic=True, use_cell_lists=True,
                               eps=self.eps, sca=self.sca,
                               radii=self.radii, boxvec=self.boxvec,
                               reference_coords=self.x_initial,
                               rcut=self.rcut, ndim=self.ndim,
                               ncellx_scale=self.ncellx_scale,
                               frozen_atoms=self.frozen_atoms2)
        self.tol = 1e-7
        self.maxstep = 1.0
        self.nstepsmax = int(1e4)
    def optimize(self, nr_samples=1):
        self.x_initial_red = []
        for a in xrange(self.N):
            if a not in self.frozen_atoms:
                self.x_initial_red.append(self.x_initial[self.ndim * a])
                self.x_initial_red.append(self.x_initial[self.ndim * a + 1])
        
        self.optimizer = ModifiedFireCPP(self.x_initial_red, self.potential, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_ = LBFGS_CPP(self.x_initial_red, self.potential)
        self.optimizer_cells = ModifiedFireCPP(self.x_initial_red, self.potential_cells, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_cells_ = LBFGS_CPP(self.x_initial_red, self.potential_cells)
        t0 = time.time()
        self.optimizer.run(self.nstepsmax)
        t1 = time.time()
        self.optimizer_cells.run(self.nstepsmax)
        t2 = time.time()
        self.optimizer_.run(self.nstepsmax)
        t3 = time.time()
        self.optimizer_cells_.run(self.nstepsmax)
        t4 = time.time()
        res_x_final = self.optimizer.get_result()
        res_x_final_cells = self.optimizer_cells.get_result()
        self.x_final = res_x_final.coords
        self.x_final_cells = res_x_final_cells.coords
        print "number of particles:", self.N
        print "time no cell lists:", t1 - t0, "sec"
        print "time cell lists:", t2 - t1, "sec"
        print "ratio:", (t1 - t0) / (t2 - t1)
        for (xci, xi) in zip(self.x_final_cells, self.x_final):
            passed = (np.abs(xci - xi) < 1e-10)
            if (passed is False):
                print "xci", xci
                print "xi", xi
                assert(passed)
        print "res_x_final.success:", res_x_final.success
        print "res_x_final_cells.success:", res_x_final_cells.success
        assert(res_x_final.success == res_x_final_cells.success)
        if not res_x_final.success:
            print "res_x_final.rms:", res_x_final.rms
            print "res_x_final.nfev:", res_x_final.nfev
        assert(res_x_final.success)
        if not res_x_final_cells.success:
            print "res_x_final_cells.rms:", res_x_final_cells.rms
            print "res_x_final_cells.nfev:", res_x_final_cells.nfev
        assert(res_x_final_cells.success)
        print "energy no cell lists:", res_x_final.energy
        print "energy cell lists:", res_x_final_cells.energy
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
    measurement(10, 10, 0.1)
    measurement_frozen(10, 10, 0.1)
