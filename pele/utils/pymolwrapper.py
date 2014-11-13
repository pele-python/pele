import pymol
from pymol import cmd, cgo
import numpy as np

def start():
    pymol.finish_launching()

def draw_spheres(coords, model, frame, radius=0.5):
    spheres=[]
    for x in coords.reshape(coords.size/3,3):
        spheres.extend([cgo.COLOR, 1.0, 0.0, 0.0])
        spheres.extend([cgo.SPHERE, x[0], x[1], x[2], radius])

    cmd.load_cgo(spheres, model, frame)

# sn402: new function for visualising rigid molecules
def draw_rigid(coords, model, frame, colour, bondslist=[], radius=0.5):
    ''' colour should be a 3-tuple representing the colour for these spheres.
        bondslist is a list of 2-tuples containing the indices of atoms which
        should be bonded together.
    '''
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
