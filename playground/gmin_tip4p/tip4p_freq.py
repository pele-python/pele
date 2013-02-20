from pygmin.thermodynamics import normalmode_frequencies
from tip4p_system import TIP4PSystem
from pygmin.angleaxis.aamindist import TransformAngleAxisCluster
from pygmin.utils import rotations
import numpy as np

system = TIP4PSystem()
db = system.create_database(db="tip4p_8.sqlite")
pot = system.get_potential()
transform = TransformAngleAxisCluster(system.aasystem)
coords = db.minima()[0].coords
#coords = db.transition_states()[0].coords
energy = pot.getEnergy(coords)
#coords = np.loadtxt("2.txt").flatten()
print energy

site = system.aasystem.sites[0]
xi = system.aasystem.sites[0].atom_positions
mi = system.aasystem.sites[0].atom_masses
M = np.sum(mi)
metric = np.identity(coords.size)

ca = system.aasystem.coords_adapter(coords)
S0 = site.Sm.copy()

gn = np.zeros([3,3])
for xx, mm in zip(xi, mi):
    y = xx # np.dot(R, xx)
    gn += mm*(np.dot(y,y)*np.identity(3) - np.outer(y,y))

# construct instantaneous metric tensor
for i, x, p in zip(xrange(coords.size/6), ca.posRigid, ca.rotRigid):
    j1 = 3*i
    j2 = 3*i + ca.nrigid*3
    metric[j1:j1+3,j1:j1+3] = np.identity(3)*M
    R = rotations.aa2mx(p)
    gx, gp = site.metric_tensor(np.zeros(3))
    #site.Sm = np.dot(R, np.dot(gp, R.transpose()))
    gp = np.dot(R, np.dot(gp, R.transpose()))
    gn = np.zeros([3,3])
    for xx, mm in zip(xi, mi):
        y = np.dot(R, xx)
        gn += mm*(np.dot(y,y)*np.identity(3) - np.outer(y,y))
        
    metric[j2:j2+3,j2:j2+3] = gn

hess = np.loadtxt('/home/vr274/tip4p_frequencies/test2inst').reshape([coords.size, coords.size])

g_inv = np.linalg.inv(metric)
A = np.dot(g_inv, hess)
print np.real(np.linalg.eigvals(A)) #[0:10] 
#print np.real(np.linalg.norm(np.dot(g_inv, metric) - np.identity(2*n)))
A = hess  #np.dot(g_inv, hess)
print "norm", np.linalg.norm(A - A.transpose())
#exit()

coords = db.minima()[0].coords

hess_read = np.loadtxt('/home/vr274/tip4p_frequencies/test1lab').reshape([coords.size, coords.size])
hess = pot.NumericalHessian(coords, eps=1e-4)

#hess = hess_read

metric = system.aasystem.metric_tensor(coords)
frq1 = normalmode_frequencies(hess, metric)

coords = np.loadtxt("2.txt").flatten()
hess_read = np.loadtxt('/home/vr274/tip4p_frequencies/test2lab').reshape([coords.size, coords.size])

hess = pot.NumericalHessian(coords, eps=1e-4)
#hess = hess_read
metric = system.aasystem.metric_tensor(coords)
frq2 = normalmode_frequencies(hess, metric)
energy = pot.getEnergy(coords)
print energy

print frq1 - frq2
A = np.dot(np.linalg.inv(metric), hess_read)
print "-------------------------"
print np.sort(np.real(frq1))
print np.sort(np.real(frq2))
#print np.real(np.linalg.eigvals(A)) #[0:10]
