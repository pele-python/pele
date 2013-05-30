import numpy as np
import exceptions as exc

class Transform(np.matrix):
    def __new__(cls, data, *args, **kwargs):
        obj = np.asarray(data).view(cls)
        return obj
    def __init__(self, *args, **kwargs):
        raise exc.NotImplementedError
    def __call__(self):
        raise exc.NotImplementedError
            
class LinearTransform(Transform):
    def __new__(cls, transform = np.identity(3), *args, **kwargs):
        obj = super(cls, LinearTransform).__new__(cls, transform)
        return obj
    def __init__(self, transform = np.identity(3), orthogonal = False, *args, **kwargs):
        self.coords_shape = (3, -1)
        self.orthogonal = orthogonal
    def __call__(self, coords):
        coords_matrix = np.reshape(coords, self.coords_shape).T
        return self * coords_matrix
    def __mul__(self, rhs):
        """
        Redefine multiplication to make the product of two orthogonal matrices orthogonal.
        We use rhs.copy(), because the implementation of np.matrix.__mul__ calls the
        asmatrix() function, re-initialising rhs (and thus setting self.orthogonality to
        False).
        """
        product = super(LinearTransform, self).__mul__(np.matrix(rhs).copy())
        if hasattr(rhs, "orthogonal"):
            product.orthogonal = (self.orthogonal & rhs.orthogonal)
        else:
            product.orthogonal = False
        return product        
    def getI(self):
        """
        This function returns the inverse matrix of LinearTransform, either by
        calculating the transpose for orthogonal matrices, or otherwise calling the 
        getI() method of the Transform (and thus Matrix) class.
        """
        if self.orthogonal == True:
            inverse = self.T
        else:
            inverse = super(LinearTransform, self).getI()
        return inverse

class AffineTransform(Transform):
    def __new__(cls,
                 linear_transform = LinearTransform(),
                 translation = np.zeros((3, 1)),
                 *args, **kwargs):
        aug_linear = np.vstack((linear_transform, np.zeros(3)))
        aug_trans  = np.vstack((np.reshape(translation, (3, 1)), np.identity(1)))
        data = np.hstack((aug_linear, aug_trans))        
        obj = super(cls, AffineTransform).__new__(cls, data)
        return obj
    def __init__(self, *args, **kwargs):
        self.coords_shape = (-1, 3)
        self.linear_transform = LinearTransform(self[0:3, 0:3])
        self.translation = [[self.flat[4 * i + 3]] for i in range(3)]
    def __call__(self, coords):
        coords_matrix = np.reshape(coords, self.coords_shape).T
        aug_coords    = np.vstack((coords_matrix, np.ones_like(coords_matrix[0])))
        return self * aug_coords
    def getI(self):
        inverse = AffineTransform(linear_transform = self.linear_transform.I,
                                  translation = -self.linear_transform.I * self.translation)
        return inverse

def identity(affine = False):
    if affine == True:
        return AffineTransform()
    else:
        return LinearTransform(orthogonal = True)

def inversion(affine = False):
    linear_inversion = LinearTransform([[-1,  0,  0],
                                        [ 0, -1,  0],
                                        [ 0,  0, -1]],
                                        orthogonal = True)
    if affine == True:
        return AffineTransform(linear_transform = linear_inversion)
    else:
        return linear_inversion
    
def proper_rotation(axis, angle, affine = False, right_handed = True):
    """
    Creates a linear or affine transformation object which corresponds to a proper
    rotation about the given axis (which passes through the origin) through the 
    given angle. By default, the rotation is right-handed (i.e. counterclockwise) 
    and returns a LinearTransform object.
    """
    # Create the matrices used in the matrix form of Rodrigues' rotation formula
    u = axis / np.linalg.norm(axis)
    identity = np.identity(3)
    tensor_prod = np.array([[u[0]*u[0], u[0]*u[1], u[0]*u[2]],
                            [u[1]*u[0], u[1]*u[1], u[1]*u[2]],
                            [u[2]*u[0], u[2]*u[1], u[2]*u[2]]])
    if right_handed == True:
        cross_prod = np.array([[   0 , -u[2],  u[1]], 
                               [ u[2],    0 , -u[0]],
                               [-u[1],  u[0],    0 ]])
    else:
        cross_prod = -1.0 * np.array([[    0, -u[2],  u[1]], 
                                      [ u[2],     0, -u[0]],
                                      [-u[1],  u[0],     0]])
    # Construct the rotation matrix from those matrices and angles
    lin_rotation_matrix = (  identity * np.cos(angle) 
                           + cross_prod * np.sin(angle)
                           + tensor_prod * (1 - np.cos(angle)))
    # Create the linear transform object
    lin_rotation_transform = LinearTransform(lin_rotation_matrix, orthogonal = True)
    if affine == False:
        transform = lin_rotation_transform
    else:
        transform = AffineTransform(linear_transform = lin_rotation_transform)
    return transform

def reflection(plane_normal, affine = False):
    """
    Creates a linear or affine transformation object which corresponds to a reflection
    in a plane with normal plane_normal and contains the origin. By default, this
    returns a LinearTransform object.
    """
    v = np.matrix(plane_normal / np.linalg.norm(plane_normal))
    identity = np.identity(3)
    reflection_matrix = identity - 2 * v.T * v
    lin_reflection_transform = LinearTransform(reflection_matrix, orthogonal = True)
    if affine == False:
        transform = lin_reflection_transform
    else:
        transform = AffineTransform(linear_transform = lin_reflection_transform)
    return transform

def improper_rotation(axis, angle, affine = False, right_handed = True):
    """
    Creates a linear or affine transformation object which corresponds to an improper
    rotation (i.e. a rotation coupled with a reflection about the axis' orthogonal
    plane) about the given axis through the given angle. This axis also passes through
    the origin. By default, the rotation is right-handed (i.e. counterclockwise) and 
    returns a LinearTransform object. 
    """
    rot = proper_rotation(axis = axis,
                          angle = angle,
                          right_handed = right_handed)
    ref = reflection(plane_normal = axis)
    lin_improper_rotation = ref * rot
    if affine == False:
        transform = lin_improper_rotation
    else:
        transform = AffineTransform(linear_transform = lin_improper_rotation)
    return transform
    

# -------------------------------------- MAIN -----------------------------------------        

if __name__ == "__main__":
#    coords = np.array([-1, 0, 1, 1, 2, 1, 1, 0, -1, -5, 3, 1])
    trans  = np.array([1, 0, 0])
#    trans2 = np.array([11, -14, 13])
    # APPLY THE TRANSLATION TO THE COORDS USING AFFINE TRANSFORM
#    test = AffineTransform(translation = trans)
#    print "Test: ", test
#    test.translation[:,0] = np.reshape(trans2, (-1,1))
#    print "Test altered:", test
#    ref = reflection(trans, affine = False)
    improper = improper_rotation(trans, np.pi * 0.5, affine = True)
    print improper.getI()
#    print test(coords)
#    print inversion(affine=False)
#    test.print_matrix()
#    test2 = LinearTransform()
#    test2.print_matrix()
#    print test2(coords)[:,2]
#    x_axis = np.array([1,0,0])
#    y_axis = np.array([0,1,0])
#    z_axis = np.array([0,0,1])
#    print proper_rotation(x_axis, np.pi / 4.0, right_handed = False)
    