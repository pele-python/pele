import numpy as np

def findrotation_kabsch(coords1, coords2):
    '''
    Kabsch, Wolfgang, (1976) "A solution of the best rotation to relate two sets of vectors", Acta Crystallographica 32:922
    '''
    
    # check if arrays are of same size
    if(coords1.size != coords2.size):
        raise BaseException("dimension of arrays does not match")
    
    # reshape the arrays
    x1 = coords1.reshape(-1,3)
    x2 = coords2.reshape(-1,3)
    
    # determine number of atoms
    natoms = x1.shape[0]
    
    # set both com to zero
    com1 = np.sum(x1,axis=0) / float(natoms)
    com2 = np.sum(x2,axis=0) / float(natoms)
    x1 -= com1
    x2 -= com2
  
    # calculate covariance matrix
    A = np.dot( x2.transpose(), x1)
    # and do single value decomposition
    u, s, v = np.linalg.svd(A)
 
    if np.linalg.det(u) * np.linalg.det(v) + 1.0 < 1e-8:
        s[-1] = -s[-1]
        u[:,-1] = -u[:,-1]
 
    return  np.dot(u, v)
    

if __name__ == "__main__":
    from pygmin.utils import rotations
    x1 = np.random.random(24)
    mx = rotations.q2mx(rotations.random_q())
    
    x2 = np.dot(mx,x1.reshape(-1,3).transpose()).transpose().reshape(-1)
    print x2-x1
    print mx-findrotation_kabsch(x1,x2)