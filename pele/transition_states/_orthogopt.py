from __future__ import print_function
import numpy as np

from ._orthogoptf import orthogopt as orthogoptf

__all__ = ["orthogopt", "orthogopt_translation_only"]


def orthogopt(v, coords, norm=False, translation_only=False):
    """
    make a vector orthogonal to eigenvectors of the Hessian corresponding to overall 
    translations and rotations.

    Parameters
    ----------
    v : numpy array
        the vector to make orthogonal
    coords : numpy array
        the zero eivengectors are calculated for this structure
    norm : bool, optional
        normalize v before returning?
    translation_only : bool
        if True, orthogonalize to translations but not rotations
    
    Returns
    -------
    v : 
        the orthogonalized vector
    """
    orthogoptf(v, coords, norm, translation_only)
    return v


def orthogopt_translation_only(v, coords, norm=False, translation_only=True):
    """a convience wrapper for orthogopt with the translation_only flag set to True"""
    return orthogopt(v, coords, norm=norm, translation_only=translation_only)


def _subcross(vec3, redcoords, n):
    """
    a utility function for orthogopt_slow
    """
    if n == 0:
        pm = 1.0
        i = 1
        j = 2
    elif n == 1:
        pm = -1.0
        i = 0
        j = 2
    else:
        pm = 1.0
        i = 0
        j = 1
    dummy1 = pm * np.sum(vec3[:, i] * redcoords[:, j] - vec3[:, j] * redcoords[:, i])
    dummy2 = np.sum((redcoords[:, i]) ** 2 + (redcoords[:, j]) ** 2)
    vdot = 0.
    if dummy2 > 0.:
        vdot = np.abs(dummy1) / np.sqrt(dummy2)
        dummy3 = dummy1 / dummy2
        # print "dummy1, dummy2", dummy1, dummy2, dummy3
        vec3[:, i] -= pm * dummy3 * redcoords[:, j]
        vec3[:, j] += pm * dummy3 * redcoords[:, i]
    return vdot


def orthogopt_slow(vec, coords, otest=False):
    """
    make vec orthogonal to eigenvectors of the Hessian corresponding to overall 
    translations and rotations.
    
    this is a really ugly recoding in python of the fortran othogopt
    """

    coords = np.reshape(coords, [-1, 3])
    natoms = len(coords[:, 0])
    com = coords.sum(0) / natoms

    vec3 = np.reshape(vec, [-1, 3])
    redcoords = coords - com

    vdot = np.zeros(3)
    vdottol = 1e-6

    for ncheck in range(100):
        vdot[:] = 0.
        ncheck += 1

        for i in range(3):
            veccom = vec3[:, i].sum() / natoms
            vdot[i] = veccom * np.sqrt(float(natoms))
            # print "vdot translation", vdot[i], i
            vec3[:, i] -= veccom
            if otest: vec3 /= np.linalg.norm(vec3)

        if False:
            if np.max(vdot) > vdottol:
                continue

        for i in range(3):
            vdot[i] = _subcross(vec3, redcoords, i)
            # print "vdot rotation   ", vdot[i], i
            if otest: vec3 /= np.linalg.norm(vec3)

        if np.max(vdot) <= vdottol:
            break

    if np.max(vdot) > vdottol:
        print("WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT")
        print("         max(vdot)", np.max(vdot))

    vec = np.reshape(vec3, [-1])
    return vec


if __name__ == "__main__":
    np.random.seed(0)
    natoms = 10
    x = np.random.rand(3 * natoms)
    v = np.random.rand(3 * natoms)
    vold = v.copy()
    v1 = orthogopt(v, x, True)
    v2 = orthogopt_slow(vold, x, True)
    print(v1 - v2)
    print("max difference between two methods", np.max(np.abs(v1 - v2)))
    

