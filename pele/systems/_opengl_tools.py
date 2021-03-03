import numpy as np

try:
    from OpenGL import GL, GLUT
except ImportError:
    pass


def subtract_com(coords):
    com = np.mean(coords, axis=0)
    return coords - com[np.newaxis, :]


def draw_sphere(xyz, radius=0.5, color=None):
    if color is not None:
        change_color(color)
    GL.glPushMatrix()
    GL.glTranslate(xyz[0], xyz[1], xyz[2])
    GLUT.glutSolidSphere(radius, 30, 30)
    GL.glPopMatrix()


def change_color(color):
    GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)


def draw_atoms(coords, atomlist, color=None, radius=0.5):
    coords = coords.reshape(-1, 3)
    if color is not None:
        change_color(color)
    for i in atomlist:
        draw_sphere(coords[i, :], radius=radius)


def draw_atomic_single_atomtype(coords, index, subtract_com=False, radius=0.5):
    """
    tell the gui how to represent your system using openGL objects
    
    Parameters
    ----------
    coords : array
    index : int
        we can have more than one molecule on the screen at one time.  index tells
        which one to draw.  They are viewed at the same time, so they should be
        visually distinct, e.g. different colors.  accepted values are 1 or 2        
    """
    coords = coords.reshape(-1, 3)
    if subtract_com:
        com = np.mean(coords, axis=0)
        coords = coords - com[np.newaxis, :]
    for x in coords:
        draw_sphere(x, radius=radius)


def draw_atomic_binary(coordslinear, index, Aatoms, Batoms, subtract_com=False,
                       rA=0.5, rB=0.44):
    """
    tell the gui how to represent your system using openGL objects
    
    Parameters
    ----------
    coords : array
    index : int
        we can have more than one molecule on the screen at one time.  index tells
        which one to draw.  They are viewed at the same time, so they should be
        visually distinct, e.g. different colors.  accepted values are 1 or 2        
    """
    coords = coordslinear.reshape(-1, 3)
    if subtract_com:
        com = np.mean(coords, axis=0)
        coords = coords - com[np.newaxis, :]

    if index == 1:
        color = [0.65, 0.0, 0.0, 1.]
    else:
        color = [0.00, 0.65, 0., 1.]
    draw_atoms(coords, Aatoms, color, radius=rA)

    if index == 1:
        color = [0.25, 0.00, 0., 1.]
    else:
        color = [0.00, 0.25, 0., 1.]

    draw_atoms(coords, Batoms, color, radius=rB)


def draw_atomic_binary_polydisperse(coordslinear, index, bdim=3, subtract_com=False, radii=None, Batoms=None):
    """
    tell the gui how to represent your system using openGL objects
    
    Parameters
    ----------
    coords : array
    bdim : box dimension
    index : int
        we can have more than one molecule on the screen at one time.  index tells
        which one to draw.  They are viewed at the same time, so they should be
        visually distinct, e.g. different colors.  accepted values are 1 or 2
    Batoms: list of atoms of type B        
    """
    assert (radii is not None)
    if Batoms is None:
        Batoms = np.ones(len(coordslinear) / bdim)

    if bdim == 2:
        # insert 0 every 2 coordinates
        j = 0
        coordslinear = coordslinear.tolist()
        for i in range(2, len(coordslinear) + 1, 2):
            coordslinear.insert(i + j, 0.0)
            j += 1
        coordslinear = np.array(coordslinear)

    coords = coordslinear.reshape(-1, 3)

    if subtract_com:
        com = np.mean(coords, axis=0)
        coords = coords - com[np.newaxis, :]

    for i, _ in enumerate(coords):
        if Batoms[i] == 1:
            if index == 1:
                color = [0.65, 0.0, 0.0, 1.]
            else:
                color = [0.00, 0.65, 0., 1.]
        else:
            if index == 1:
                color = [0.25, 0.00, 1., 1.]
            else:
                color = [0.00, 0.25, 1., 1.]

        draw_atoms(coords, [i], color, radius=radii[i])


def draw_cone(X1, X2, rbase=0.1, rtop=0.0, color=None):
    """draw a cylinder from X1 to X2"""
    if color is not None:
        change_color(color)
    from OpenGL import GL, GLU

    z = np.array([0., 0., 1.])  # default cylinder orientation
    p = X2 - X1  # desired cylinder orientation
    r = np.linalg.norm(p)
    t = np.cross(z, p)  # angle about which to rotate
    a = np.arccos(np.dot(z, p) / r)  # rotation angle
    a *= (180. / np.pi)  # change units to angles
    GL.glPushMatrix()
    GL.glTranslate(X1[0], X1[1], X1[2])
    GL.glRotate(a, t[0], t[1], t[2])
    g = GLU.gluNewQuadric()
    GLU.gluCylinder(g, rbase, rtop, r, 30, 30)  # I can't seem to draw a cylinder
    GL.glPopMatrix()


def draw_cylinder(X1, X2, radius=.1, color=None):
    draw_cone(X1, X2, rbase=radius, rtop=radius, color=color)

def draw_box(boxvec, radius=0.05):
    """draw the edges of a box with center at the origin"""
    e = np.eye(3)
    e[0,0] = boxvec[0]
    e[1,1] = boxvec[1]
    e[2,2] = boxvec[2]
    
    from itertools import product
    corners = [np.array(x) for x in product([0,1], repeat=3)]
    
    x0 = - boxvec / 2
    for i, c1 in enumerate(corners):
        for c2 in corners[:i]:
            if np.sum(np.abs(c1-c2)) == 1:
                # these corners are adjacent, i.e not diagonal
                x1 = x0 + np.dot(e, c1)
                x2 = x0 + np.dot(e, c2)
                draw_cylinder(x1, x2, radius=radius)

    

