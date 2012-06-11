import pymol
from pymol import cmd, cgo

def start():
    pymol.finish_launching()

def draw_spheres(coords, model, frame, radius=0.5):
    spheres=[]
    for x in coords.reshape(coords.size/3,3):
        spheres.extend([cgo.COLOR, 1.0, 0.0, 0.0])
        spheres.extend([cgo.SPHERE, x[0], x[1], x[2], radius])

    cmd.load_cgo(spheres, model, frame)
