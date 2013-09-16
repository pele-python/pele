import numpy as np
try:
    from OpenGL import GL, GLUT
except ImportError:
    pass

def subtract_com(coords):
    com = np.mean(coords, axis=0)
    return coords - com[np.newaxis,:]

def draw_sphere(xyz, radius=0.5):
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
        draw_sphere(coords[i,:])

def draw_atomic_single_atomtype(coords, index, subtract_com=False):
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
    from OpenGL import GL,GLUT
    coords = coords.reshape(-1, 3)
    if subtract_com:
        com = np.mean(coords, axis=0)
        coords = coords - com[np.newaxis,:]
    for x in coords:
        draw_sphere(x)


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
    from OpenGL import GL,GLUT
    coords = coordslinear.reshape(-1, 3)
    if subtract_com:
        com = np.mean(coords, axis=0)
        coords = coords - com[np.newaxis,:]
                          
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
