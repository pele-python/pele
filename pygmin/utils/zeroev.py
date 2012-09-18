'''
Created on 2 Aug 2012

@author: ruehle
'''

import numpy as np

def for_cluster(coords):
        # translational eigenvectors
        x1 = np.zeros(coords.shape)
        x2 = x1.copy()
        x3 = x1.copy()
        x1.reshape(coords.size/3,3)[:,0] = 1.
        x2.reshape(coords.size/3,3)[:,1] = 1.
        x3.reshape(coords.size/3,3)[:,2] = 1.
        Rx = np.array([[ 0.,  0.,  0.],
                       [ 0.,  0.,  1.],
                       [ 0., -1.,  0.]])
        Ry = np.array([[ 0.,  0.,  1.],
                       [ 0.,  0.,  0.],
                       [-1.,  0.,  0.]])
        Rz = np.array([[ 0.,  1., 0.],
                       [-1.,  0., 0.],
                       [ 0.,  0., 0.]])
        x = coords.reshape(coords.size/3,3)
        com = x.sum(0)/x.shape[0]
        #print com
        r1 = np.dot(Rx,(x-com).transpose()).transpose().reshape(coords.shape)
        r2 = np.dot(Ry,(x-com).transpose()).transpose().reshape(coords.shape)
        r3 = np.dot(Rz,(x-com).transpose()).transpose().reshape(coords.shape)
        
        #for u,v in zip(r2.reshape(coords.size/3,3), x.reshape(coords.size/3,3)):
        #    print np.dot(u,v-com)
        return [x1/np.linalg.norm(x1), x2/np.linalg.norm(x2), x3/np.linalg.norm(x3), 
                r1/np.linalg.norm(r1), r2/np.linalg.norm(r2), r3/np.linalg.norm(r3)]
        
def gramm_schmidt(v):
    u = []
    for vi in v:
        vn=vi.copy()
        for uk in u:
            vn -= np.dot(vn,uk)*uk
        #print "lz",len(u)
        vn = vi / np.linalg.norm(vi)
        u.append(vn)
    return u

if __name__ == '__main__':
    from pygmin.potentials import lj
    pot = lj.LJ()
    x = np.array([-0.,0.,0.,1.,1.,1.])
    
    #x = np.random.random(6)
    print x
    u = gramm_schmidt(for_cluster(x))
    for i in u:
        print (pot.getEnergy(x + 1e-6) - pot.getEnergy(x))/1e-6,i
    print np.dot(u[3],u[4]),np.dot(u[3],u[5]),np.dot(u[5],u[4])
    print u[5],u[4]