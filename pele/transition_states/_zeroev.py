"""
Created on 2 Aug 2012

@author: ruehle
"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as np

__all__ = ["zeroEV_translation", "zeroEV_rotation", "zeroEV_cluster", "gramm_schmidt"]


def zeroEV_translation(coords):
    """
    return the set of zero eigenvalue vectors corresponding to uniform translation
    """
    x1 = np.zeros(coords.shape)
    x2 = np.zeros(coords.shape)
    x3 = np.zeros(coords.shape)
    x1.reshape(-1, 3)[:, 0] = 1.
    x2.reshape(-1, 3)[:, 1] = 1.
    x3.reshape(-1, 3)[:, 2] = 1.
    return [x1 / np.linalg.norm(x1), x2 / np.linalg.norm(x2), x3 / np.linalg.norm(x3)]


def zeroEV_rotation(coords):
    """
    return the set of zero eigenvalue vectors corresponding to global rotation
    
    note that these are not necessarily orthogonal.  use `gramm_schmidt` to make them orthogonal
    
    See Also
    --------
    gramm_schmidt
    """
    # translational eigenvectors
    Rx = np.array([[0., 0., 0.],
                   [0., 0., 1.],
                   [0., -1., 0.]])
    Ry = np.array([[0., 0., 1.],
                   [0., 0., 0.],
                   [-1., 0., 0.]])
    Rz = np.array([[0., 1., 0.],
                   [-1., 0., 0.],
                   [0., 0., 0.]])
    x = coords.reshape(coords.size // 3, 3)
    com = x.sum(0) / x.shape[0]

    r1 = np.dot(Rx, (x - com).transpose()).transpose().reshape(coords.shape)
    r2 = np.dot(Ry, (x - com).transpose()).transpose().reshape(coords.shape)
    r3 = np.dot(Rz, (x - com).transpose()).transpose().reshape(coords.shape)

    return [r1 / np.linalg.norm(r1), r2 / np.linalg.norm(r2), r3 / np.linalg.norm(r3)]


def zeroEV_cluster(coords):
    """
    return the set of zero eigenvalue vectors corresponding to uniform translation and rotation
    """
    return zeroEV_translation(coords) + zeroEV_rotation(coords)


def gramm_schmidt(v):
    """
    make a set of vectors orthogonal to each other.
    """
    u = []
    for vi in v:
        vn = vi.copy()
        for uk in u:
            vn -= np.dot(vn, uk) * uk
        vn /= np.linalg.norm(vn)
        u.append(vn)
    return u


def orthogonalize(v, ozev):
    """make vector v orthogonal to all vectors in list ozev
    
    Parameters
    ----------
    v : array
        vector to orthogonalize
    ozev : list of arrays
        v will be made orthogonal to all of these vectors. 
    """
    for u in ozev:
        v -= np.dot(v, u) * u
    return v


#
# testing only below here
#

def test():  # pragma: no cover
    from ._orthogopt import orthogopt_slow, orthogopt

    natoms = 105
    for i in range(1):
        x = np.random.random(3 * natoms) * 5
        xx = x.reshape(-1, 3)
        com = xx.sum(0) / xx.shape[0]
        xx -= com
        v = np.random.random(3 * natoms)
        test1 = orthogopt_slow(v.copy(), x.copy())
        ozev = gramm_schmidt(zeroEV_cluster(x))

        orthogonalize(v, ozev)

        print(np.linalg.norm(v - test1))
    exit()

    from pele.potentials import lj

    pot = lj.LJ()
    x = np.array([-1., 0., 0., 1., 0., 0., 0., 1., 1., 0., -1., -1.])
    x = np.random.random(x.shape)
    print(x)
    v = zeroEV_cluster(x)
    print(np.dot(v[0], v[1]), np.dot(v[0], v[2]), np.dot(v[1], v[2]))
    print(np.dot(v[3], v[4]), np.dot(v[3], v[5]), np.dot(v[5], v[4]))
    u = gramm_schmidt(zeroEV_cluster(x))
    for i in u:
        print((pot.getEnergy(x + 1e-4 * i) - pot.getEnergy(x)) / 1e-4, i)
    print(np.dot(u[3], u[4]), np.dot(u[3], u[5]), np.dot(u[5], u[4]))
    print("########################")

    r = np.random.random(x.shape)
    print(orthogopt(r.copy(), x.copy()) - orthogonalize(r.copy(), u))


if __name__ == '__main__':
    test()

