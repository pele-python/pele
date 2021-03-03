from . import pymol
from .pymol import cmd, cgo
import numpy as np

def start():
    pymol.finish_launching()

def draw_spheres(coords, model, frame, radius=0.5):
    spheres=[]
    for x in coords.reshape(coords.size/3,3):
        spheres.extend([cgo.COLOR, 1.0, 0.0, 0.0])
        spheres.extend([cgo.SPHERE, x[0], x[1], x[2], radius])

    cmd.load_cgo(spheres, model, frame)


def draw_rigid(coords, model, frame, colour, bondslist=[], radius=0.5):
    """ 
    Use pymol to draw a system of rigid body fragments
    
    Parameters
    ----------
    
    colour: 3-tuple 
        RBG colour for the spheres being drawn.
    bondslist: list of 2-tuples, optional
        List of atom pairs between which bonds should be drawn.
    """
    spheres = []
   
    
    for x in coords.reshape(coords.size/3,3):
        spheres.extend([cgo.COLOR, colour[0], colour[1], colour[2]])
        spheres.extend([cgo.SPHERE, x[0], x[1], x[2], radius])           

    coords = coords.reshape(coords.size/3,3)
    Rcyl = .1     
    
    for i in bondslist:
        spheres.extend([cgo.CYLINDER, coords[i[0]][0], coords[i[0]][1], coords[i[0]][2],
                        coords[i[1]][0], coords[i[1]][1], coords[i[1]][2],
                        Rcyl , 255., 255., 255. , 0., 0., 0.])
        
       
    

    cmd.load_cgo(spheres, model, frame)  
    
def draw_box(boxvec, model, frame):
    """
    Draw a box around the system to easily visualise periodic boundaries
    
    Parameters
    ----------
    boxvec: np.array
        The dimensions of the periodic box
    """  
    box = []
    Rcyl = .1    
    
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], -0.5*boxvec[1], -0.5*boxvec[2], 0.5*boxvec[0], -0.5*boxvec[1], 
                -0.5*boxvec[2], Rcyl, 1., 0., 0., 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], -0.5*boxvec[1], -0.5*boxvec[2], -0.5*boxvec[0], 0.5*boxvec[1], 
                -0.5*boxvec[2], Rcyl , 0., 1., 0. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], -0.5*boxvec[1], -0.5*boxvec[2], -0.5*boxvec[0], -0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 0., 0., 1. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, 0.5*boxvec[0], -0.5*boxvec[1], -0.5*boxvec[2], 0.5*boxvec[0], 0.5*boxvec[1], 
                -0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, 0.5*boxvec[0], -0.5*boxvec[1], -0.5*boxvec[2], 0.5*boxvec[0], -0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], 0.5*boxvec[1], -0.5*boxvec[2], 0.5*boxvec[0], 0.5*boxvec[1], 
                -0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], 0.5*boxvec[1], -0.5*boxvec[2], -0.5*boxvec[0], 0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], -0.5*boxvec[1], 0.5*boxvec[2], 0.5*boxvec[0], -0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, -0.5*boxvec[0], -0.5*boxvec[1], 0.5*boxvec[2], -0.5*boxvec[0], 0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, 0.5*boxvec[0], 0.5*boxvec[1], 0.5*boxvec[2], 0.5*boxvec[0], 0.5*boxvec[1], 
                -0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, 0.5*boxvec[0], 0.5*boxvec[1], 0.5*boxvec[2], -0.5*boxvec[0], 0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])
    box.extend([cgo.CYLINDER, 0.5*boxvec[0], 0.5*boxvec[1], 0.5*boxvec[2], 0.5*boxvec[0], -0.5*boxvec[1], 
                0.5*boxvec[2], Rcyl , 255., 255., 255. , 0., 0., 0.])                                              

    cmd.load_cgo(box, model, frame)  

