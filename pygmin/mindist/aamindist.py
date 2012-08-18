import numpy as np

def aareduce( v1 ):
    """
    make the magnitude of angle-axis vector v1 as close to zero as possible by
    applying symmetry operations
    """
    v1n = np.linalg.norm( v1 )
    n = np.round(v1n / (2.*np.pi))
    v1[:] *= (v1n-n*2.*np.pi) / v1n

def aadistance( v1, v2 ):
    """
    minimize the "distance" between two angle axis vectors
    
    perform symmetry operations on angle axis vectors v1 and v2 such that the
    cartesian distance norm(v1(:) - v2(:)) between them is minimized.
     
    the possible symmetry operations are

    v1 -> v1 * ( 1 + n * 2*pi / norm(v1) ) 

    where n is an integer
    
    modify v1 and v2 in place
    """
    #v1 = np.copy(v1_in)
    #v2 = np.copy(v2_in)

    #dold = sqrt(sum( (v2 - v1)**2 )) !for debugging

    #make the magnitude of v1 as close to zero as possible
    aareduce( v1 )
    aareduce( v2 )

    """
     We still need to check the two possibilities

     v1 -> v1 * (v1n - 2*pi)/v1n   
     v2 -> v2

     and

     v1 -> v1
     v2 -> v2 * (v2n - 2*pi)/v2n

     Other symmetry operations can be ruled out as options to give to closer
     vectors.  We can simplify the calculation by noticing that the smaller of
     v1 and v2 will always remain unchanged.  Finally, we can check if the
     larger vector needs to change by checking the condition (asuming v1 is the
     larger)

     norm(v1) - dot_product(v1, v2) / norm(v1)  > Pi
    """

    #dnew = sqrt(sum( (v2 - v1)**2 )) ! for debugging

    v1n = np.linalg.norm(v1)
    v2n = np.linalg.norm(v2)
    if v1n > v2n:
        d = v1n - np.dot(v1, v2) / v1n
        if d > np.pi:
            v1 *= (1.0 - 2.*np.pi / v1n)
    else:
        d = v2n - np.dot(v1, v2) / v2n
        if d > np.pi:
            v2 *= (1.0 - 2.*np.pi / v2n)
    #dnewnew = sqrt(sum( (v2 - v1)**2 )) !js850>
    #if (dold > 2) write(*,*) "js850>", dold, dnew, dnewnew !js850>
    #if (dnewnew > 3) write(*,*) "js850> something may have gone wrotg"


    return v1, v2


import unittest
class aaDistTest(unittest.TestCase):
    def setUp(self):
        self.v1 = np.random.uniform(-8,8.,[3])
        self.v2 = np.random.uniform(-8,8.,[3])
        self.v1in = np.copy(self.v1)
        self.v2in = np.copy(self.v2)
        
    def test_dist(self):
        aadistance(self.v1, self.v2)
        
        v1n = np.linalg.norm(self.v1)
        v2n = np.linalg.norm(self.v2)
        
        v1inn = np.linalg.norm(self.v1in)
        v1inn = np.linalg.norm(self.v1in)

        distbefore = np.linalg.norm(self.v1in - self.v2in)
        distafter = np.linalg.norm(self.v1 - self.v2)

        self.assertTrue(distbefore >= distafter, "distance between aa vectors increased")
        
        import pygmin.utils.rotations as rot
        
        rot1 = rot.aa2mx( self.v1 )
        rot1i = rot.aa2mx( self.v1in )
        self.assertTrue( all( (rot1 == rot1i)[:].tolist() ), "aa vector returned different rotation matrix" )
        
        rot2 = rot.aa2mx( self.v2 )
        rot2i = rot.aa2mx( self.v2in )
        self.assertTrue( all( (rot2 == rot2i)[:].tolist() ), "aa vector returned different rotation matrix" )



        

        

def test(v1in, v2in, aadistfunc = aadistance):
    v1 = np.copy(v1in)
    v2 = np.copy(v2in)

    print v1, v2
    v1n = np.linalg.norm(v1)
    v2n = np.linalg.norm(v2)
    print "initial norms", v1n, v2n
    
    aadistfunc(v1, v2)
    

    
    v3n = np.linalg.norm(v1)
    v4n = np.linalg.norm(v2)
    print "final norms  ", v3n, v4n
    import pygmin.utils.rotations as rot
    print "v1 unchanged:", all( (rot.aa2mx( v1) ==  rot.aa2mx( v1in))[:].tolist() )
    print "v2 unchanged:", all( (rot.aa2mx( v2) ==  rot.aa2mx( v2in))[:].tolist() )

    print "cartesian distance between vectors before", np.linalg.norm(v1in-v2in)
    print "cartesian distance between vectors after ", np.linalg.norm(v1-v2)

    


if __name__ == "__main__":
    v1 = np.random.uniform(-8,8.,[3])
    v2 = np.random.uniform(-8,8.,[3])

    test(v1, v2)
    
    try:
        from aadistance_fort import aadistance as aadistfort
        print ""
        print ""
        print "comparing with fortran version"
        test(v1, v2, aadistfort)
    except ImportError:
        pass
    
    unittest.main()
