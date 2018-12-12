from __future__ import absolute_import
import tempfile

from pele.systems import AtomicCluster
from pele.potentials import LJ
from pele.utils.xyz import write_xyz

__all__ = ["LJCluster"]


class LJCluster(AtomicCluster):
    """
    define the System class for a Lennard-Jones cluster
    
    Parameters
    ----------
    natoms : int 
    
    See Also
    --------
    BaseSystem, AtomicCluster
    """

    def __init__(self, natoms):
        super(LJCluster, self).__init__()
        self.natoms = natoms

        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0

    # self.params.double_ended_connect.NEBparams.reinterpolate = 1

    def get_permlist(self):
        return [list(range(self.natoms))]

    def get_potential(self):
        return LJ()

    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    potential="LJ cluster",
        )

    #
    # below here is stuff only for the gui
    #

    def draw(self, coordslinear, index, subtract_com=True):  # pragma: no cover
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
        from ._opengl_tools import draw_atomic_single_atomtype
        draw_atomic_single_atomtype(coordslinear, index, subtract_com=subtract_com)

    def load_coords_pymol(self, coordslist, oname, index=1):  # pragma: no cover
        """load the coords into pymol
        
        the new object must be named oname so we can manipulate it later
                        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object must be named oname so it can be manipulated
            later
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2
        
        Notes
        -----
        the implementation here is a bit hacky.  we create a temporary xyz file from coords
        and load the molecule in pymol from this file.  
        """
        # pymol is imported here so you can do, e.g. basinhopping without installing pymol
        from . import pymol

        # create the temporary file
        suffix = ".xyz"
        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        fname = f.name

        # write the coords into the xyz file
        from pele.mindist import CoMToOrigin

        for coords in coordslist:
            coords = CoMToOrigin(coords.copy())
            write_xyz(f, coords, title=oname, atomtypes=["LA"])
        f.flush()

        # load the molecule from the temporary file
        pymol.cmd.load(fname)

        # get name of the object just create and change it to oname
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        pymol.cmd.set_name(objectname, oname)

        # set the representation
        pymol.cmd.hide("everything", oname)
        pymol.cmd.show("spheres", oname)

        # set the color according to index
        if index == 1:
            pymol.cmd.color("red", oname)
        else:
            pymol.cmd.color("gray", oname)


#
# only for testing below here
#

def run():  # pragma: no cover
    # create the system object
    sys = LJCluster(15)

    # create a database
    db = sys.create_database()

    # do a short basinhopping run
    bh = sys.get_basinhopping(database=db, outstream=None)
    while len(db.minima()) < 2:
        bh.run(100)

    # try to connect the lowest two minima
    min1, min2 = db.minima()[:2]
    connect = sys.get_double_ended_connect(min1, min2, db)
    connect.connect()


if __name__ == "__main__":
    run()

