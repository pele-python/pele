"""
.. currentmodule:: pele.potentials

Potentials (`pele.potentials`)
=================================

This module exports native potential bindings and Python implementations.
The reference below groups the available models and wrappers by purpose.

Potential interface
-------------------

Python potentials can derive from these base classes.

.. autosummary::
   :toctree: generated/

    BasePotential
    ~_pele.PairwisePotentialInterface
    BasePotentialAtomistic

When creating a potential, implement ``getEnergy()``. Many routines also use
gradients, so implement ``getEnergyGradient()`` to avoid calculating them
numerically. See :doc:`/tutorial_potential` for a worked example.

Pair potentials
---------------

.. autosummary::
    :toctree: generated/

    LJ
    BLJCut
    LJCut
    LJpshift
    Frenkel
    Morse
    HS_WCA
    WCA
    WCANeighborList
    WCAAtomList

Inverse-power potentials
------------------------

.. autosummary::
    :toctree: generated/

    InversePower
    InversePowerHS
    InversePowerStillinger
    InversePowerStillingerCut
    InversePowerStillingerCutQuad
    PyInversePower

Many-body potentials
--------------------

.. autosummary::
    :toctree: generated/

    ~ATLJ.ATLJ
    ATLJCPP

Harmonic and Gaussian potentials
--------------------------------

.. autosummary::
    :toctree: generated/

    Harmonic
    RadialGaussian
    SumGaussianPot

Spin potentials
---------------

.. autosummary::
    :toctree: generated/

    XYModel
    HeisenbergModel
    HeisenbergModelRA
    MeanFieldPSpinSpherical

Potential wrappers
------------------

.. autosummary::
    :toctree: generated/

    CombinedPotential
    FrozenPotentialWrapper

Cost and test functions
-----------------------

.. autosummary::
    :toctree: generated/

    MLCost
    PoweredCosineSum

External potentials
-------------------

.. autosummary::
   :toctree: generated/

    GMINPotential

``LAMMPSPotential`` is available when the optional LAMMPS dependency is
installed. See the `Python adapter source
<https://github.com/pele-python/pele/blob/main/pele/potentials/lammps_potential_python.py>`_
and `native adapter source
<https://github.com/pele-python/pele/blob/main/pele/potentials/lammps_pele_cython.pyx>`_
for its energy and gradient interface.

"""



from .potential import *
from ._frozen_dof import FrozenPotentialWrapper

# from lj import *
from ._lj_cpp import LJ, BLJCut
from ._frenkel import Frenkel
from ._hs_wca_cpp import HS_WCA
from ._inversepower_cpp import InversePower
from ._inversepower_stillinger_cpp import InversePowerStillinger
from ._inversepower_stillinger_cut_cpp import InversePowerStillingerCut
from ._inversepower_stillinger_cut_quad import InversePowerStillingerCutQuad
from .combine_potentials import CombinedPotential
from ._inversepower_hs_cpp import InversePowerHS
from ._wca_cpp import *
from ._harmonic_cpp import Harmonic
from ._radial_gaussian_cpp import RadialGaussian
from ._sumgaussianpot_cpp import SumGaussianPot
from ._pspin_spherical_cpp import MeanFieldPSpinSpherical
from .ATLJ import *
from .atlj import ATLJ as ATLJCPP
from .gminpotential import *
from .heisenberg_spin import *
from .heisenberg_spin_RA import *
from .ljpshiftfast import *
from .ljcut import *

# from potential import *
# from salt import *
# from soft_sphere import *
# from stockmeyer import *
from .xyspin import *
from .morse import Morse
from .ml import MLCost
from .inverse_power_py import PyInversePower
from .cpp_test_functions import PoweredCosineSum

try:
    import lammps
    have_lammps = True
    try:
        from pele.potentials.lammps_pele_cython import LAMMPSPotential
    except ImportError:
        from pele.potentials.lammps_potential_python import LAMMPSPotential
except ImportError:
    have_lammps = False
