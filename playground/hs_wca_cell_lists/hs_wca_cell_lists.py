from __future__ import division
import numpy as np
from pele.potentials import HS_WCA, HS_WCAPeriodicCellLists
from pele.optimize import ModifiedFireCPP, LBFGS_CPP
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

def save_pdf(plt, file_name):
    pdf = PdfPages(file_name)
    plt.savefig(pdf, format="pdf")
    pdf.close()
    plt.close()
    
def to_string(input, digits_after_point = 16):
    format_string = "{0:."
    format_string += str(digits_after_point)
    format_string += "f}"
    return format_string.format(input)

class Config2D(object):
    def __init__(self, nparticles_x, amplitude):
        self.LX = nparticles_x
        self.LY = self.LX
        self.nparticles_x = nparticles_x
        self.N = self.nparticles_x ** 2
        self.amplitude = amplitude
        self.x = np.zeros(2 * self.N)
        for particle in xrange(self.N):
            pid = 2 * particle
            self.x[pid] = particle % self.LX
            self.x[pid + 1] = int(particle / self.LX)
        self.x_initial = [xi + np.random.uniform(- self.amplitude, self.amplitude) for xi in self.x]
        self.radius = 0.25
        self.sca = 1.1
        self.radii = np.ones(self.N) * self.radius
        self.eps = 1
        self.boxvec = np.array([self.LX, self.LY])
        self.potential = HS_WCA(self.eps, self.sca, self.radii, ndim=2, boxvec=self.boxvec)
        self.rcut = 2 * (1 + self.sca) * self.radius
        self.ncellx_scale = 1
        self.potential_cells = HS_WCAPeriodicCellLists(self.eps, self.sca, self.radii, self.boxvec, self.x_initial, self.rcut, ndim = 2, ncellx_scale = self.ncellx_scale)
        self.tol = 1e-7
        self.maxstep = 1
        self.nstepsmax = 1e5
    def optimize(self, nr_samples = 1):
        self.optimizer = ModifiedFireCPP(self.x_initial, self.potential, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_ = LBFGS_CPP(self.x_initial, self.potential)
        self.optimizer_cells = ModifiedFireCPP(self.x_initial, self.potential_cells, tol = self.tol, maxstep = self.maxstep)
        self.optimizer_cells_ = LBFGS_CPP(self.x_initial, self.potential_cells)
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
        assert(res_x_final.success)
        assert(res_x_final_cells.success)
        for (xci, xi) in zip(self.x_final_cells, self.x_final):
            assert(np.abs(xci - xi) < 1e-10 * self.N)
        print "energy no cell lists:", res_x_final.energy
        print "energy cell lists:", res_x_final_cells.energy
        self.t_ratio = (t1 - t0) / (t2 - t1)
        self.t_ratio_lbfgs = (t3 - t2) / (t4 - t3)

def ll(x, a, b):
    return a * x ** b

if __name__ == "__main__":
    nr_samples = 10
    LXmax = 19
    nr_particles = []
    time_ratio = []
    time_ratio_lbfgs = []
    for LX in xrange(3, LXmax + 1):
        print "---LX: ", LX
        t_fire = 0
        t_lbfgs = 0
        nrp = 0
        for tmp in xrange(nr_samples):
            conf = Config2D(LX, 0.03)
            conf.optimize()
            t_fire = t_fire + conf.t_ratio / nr_samples
            t_lbfgs = t_lbfgs + conf.t_ratio_lbfgs / nr_samples
            nrp = conf.N
        nr_particles.append(nrp)
        time_ratio.append(t_fire)
        time_ratio_lbfgs.append(t_lbfgs)    
        print "---done:", LX, "max: ", LXmax  
    popt, pcov = curve_fit(ll, nr_particles, time_ratio)
    plt.plot(nr_particles, time_ratio, "s", label = r"M-FIRE")
    plt.plot(nr_particles, time_ratio_lbfgs, "o", label = r"LBFGS")
    xr = np.linspace(nr_particles[0], nr_particles[-1])
    plt.plot(xr, [ll(xi, popt[0], popt[1]) for xi in xr], "k", label=r"Exponent: " + to_string(popt[1], 2))
    plt.legend(loc = 2)
    plt.xlabel(r"Number of particles")
    plt.ylabel(r"Time no cell lists / time cell lists")
    save_pdf(plt, "hs_wca_cell_lists.pdf")
