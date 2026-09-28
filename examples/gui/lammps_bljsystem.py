import sys

from PyQt5 import QtWidgets

from pele.gui import run_gui
from pele.systems import AtomicClusterLAMMPS


def get_natoms(atom_type: str) -> int:
    dialog = QtWidgets.QInputDialog()
    dialog.setLabelText(f"number of {atom_type} atoms")
    dialog.setWindowTitle("Create new *LAMMPS* Lennard-Jones system")
    dialog.setInputMode(1)
    dialog.setIntMinimum(1)
    dialog.exec_()
    if not dialog.result():
        sys.exit()
    return dialog.intValue()


def lammps_blj_example():
    from lammps import lammps

    app = QtWidgets.QApplication(sys.argv)
    natoms_a = get_natoms("A")
    natoms_b = get_natoms("B")

    lmp = lammps(cmdargs=f"-screen none -log none -v na {natoms_a} -v nb {natoms_b}".split())
    lmp.file("./lammps_bljsystem.in")
    system = AtomicClusterLAMMPS(lmp)
    run_gui(system, application=app)


if __name__ == "__main__":
    lammps_blj_example()
