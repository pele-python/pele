import sys

from PyQt5 import QtWidgets

from pele.gui import run_gui
from pele.systems import AtomicClusterLAMMPS


def lammps_lj_example():
    from lammps import lammps

    app = QtWidgets.QApplication(sys.argv)

    dialog = QtWidgets.QInputDialog()
    dialog.setLabelText("number of atoms")
    dialog.setWindowTitle("Create new *LAMMPS* Lennard-Jones system")
    dialog.setInputMode(1)
    dialog.setIntMinimum(2)
    dialog.setIntValue(13)
    dialog.exec_()
    if not dialog.result():
        sys.exit()
    natoms = dialog.intValue()

    lmp = lammps(cmdargs=f"-screen none -log none -v n {natoms}".split())
    lmp.file("./lammps_ljsystem.in")
    system = AtomicClusterLAMMPS(lmp)
    run_gui(system, application=app)


if __name__ == "__main__":
    lammps_lj_example()
