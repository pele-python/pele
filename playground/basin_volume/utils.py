from __future__ import division
import numpy as np
import os
from scipy.special import gamma

def volume_nball(radius,n):
    volume = 2*np.power(np.pi,n/2)*np.power(radius,n)/(n*gamma(n/2))
    return volume

def cround(r):
    if r > 0.0:
        r = np.floor(r + 0.5)
    else:
        r = np.ceil(r - 0.5)
    return r

def trymakedir(path):
    """this function deals with common race conditions"""
    while True:
        if not os.path.exists(path): 
            try:
                os.makedirs(path)
                break
            except OSError, e:
                if e.errno != 17:
                    raise
                # time.sleep might help here
                pass
        else:
            break

def put_in_box(x, boxvec):
    x = x.reshape(-1, len(boxvec))
    x -= boxvec * np.round(x / boxvec)
    
def read_xyzd(fname):
    coords = []
    radii = []
    f = open(fname, "r")
    while True:
        xyzd = f.readline()
        if not xyzd: break
        x, y, z, d = xyzd.split()
        coords.extend([float(x),float(y),float(z)])
        radii.extend([float(d)/2])
    return np.array(coords), np.array(radii)

def read_xyzdr(fname, bdim=3):
    coords = []
    radii = []
    rattlers = []
    f = open(fname, "r")
    while True:
        xyzdr = f.readline()
        if not xyzdr: break
        x, y, z, d, r = xyzdr.split()
        coords.extend([float(x),float(y),float(z)])
        radii.extend([float(d)/2])
        for _ in xrange(bdim): 
            rattlers.extend([float(r)])
    return np.array(coords), np.array(radii), np.array(rattlers)
