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
    coords_shape = (-1, 3)
    def __new__(cls, transform = np.identity(3), *args, **kwargs):
        obj = super(cls, LinearTransform).__new__(cls, transform)
        return obj
    def __init__(self, transform = np.identity(3), orthogonal = False, *args, **kwargs):
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
    coords_shape = (-1, 3)
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
        self.linear_transform = LinearTransform(self[0:3, 0:3])
        self.translation = [[self.flat[4 * i + 3]] for i in range(3)]
    def __call__(self, coords):
        coords_matrix = np.reshape(coords, self.coords_shape).T
        aug_coords    = np.vstack((coords_matrix, np.ones_like(coords_matrix[0])))
        return np.reshape(np.array(((self * aug_coords)[0:3]).T), -1)
    def getI(self):
        inverse = AffineTransform(linear_transform = self.linear_transform.I,
                                  translation = -self.linear_transform.I * self.translation)
        return inverse

def identity(affine = False):
    """
    Returns the identity operation as a LinearTransform or AffineTransform object.  
    """
    if affine == True:
        return AffineTransform()
    else:
        return LinearTransform(orthogonal = True)

def inversion(affine = False):
    """
    Returns the inversion operation as a LinearTransform or AffineTransform object.
    """
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
    Returns a linear or affine transformation object which corresponds to a proper
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
        cross_prod = np.array([[   0 ,   u[2], -u[1]], 
                               [ -u[2],    0 ,  u[0]],
                               [  u[1], -u[0],    0 ]])
    else:
        cross_prod = -1.0 * np.array([[    0, u[2],  -u[1]], 
                                      [ -u[2],     0, u[0]],
                                      [u[1],  -u[0],     0]])
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
    Returns a linear or affine transformation object which corresponds to a reflection
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
    Returns a linear or affine transformation object which corresponds to an improper
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

def translation(translation):
    """
    Returns an affine transformation object which corresponds to a translation.
    """
    trans = AffineTransform(translation = np.array(translation))
    return trans    

# -------------------------------------- MAIN -----------------------------------------        

if __name__ == "__main__":
    x_axis = np.array([1,0,0])
    y_axis = np.array([0,1,0])
    z_axis = np.array([0,0,1])
    
    coords = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1])
    
    trans  = np.array([0, 4, 0])
    rot    = np.pi
    
    proper = proper_rotation(x_axis, rot, affine = True)
    improper = improper_rotation(x_axis, rot, affine = True)
    inverse = inversion(affine = True)
    trans_forward = translation(trans)
    trans_back = translation(-trans)
    dihedral_move = trans_back * proper * trans_forward
    print inverse(coords)
    print improper(coords)