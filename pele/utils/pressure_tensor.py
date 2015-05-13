"""
    Computes the pressure tensor and returns the pressure.
    The negative of the pressure tensor is often called the stress
    tensor, see Allen, Tidsley, "Computer Simulation of Liquids"
    pp. 60-61.
"""

##        double get_pressure_tensor(Array[double] &x, Array[double] &ptensor, double volume) except +

def getPressureTensor(pot, np.ndarray[double, ndim=1] x not None, volume, ndim):
    cdef np.ndarray[double, ndim=1] ptensor = np.zeros(ndim * ndim)
    p = self.thisptr.get().get_pressure_tensor(array_wrap_np(x), array_wrap_np(ptensor), volume)
    return p, ptensor
