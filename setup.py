from numpy.distutils.core import setup
from numpy.distutils.core import Extension

setup(name='PyGMIN', 
        version='0.1', 
        description="Pythin implementation of GMIN, OPTIM, and PATHSAMPLE",
        url='http://github.com/js850/minGMIN.git',
        packages=["pygmin", 
            "pygmin.potentials",
            "pygmin.gui",
            "pygmin.mindist",
            "pygmin.optimize",
            "pygmin.optimize.transition_state",
            "pygmin.printing",
            "pygmin.takestep",
            "pygmin.utils",
            "pygmin.wham"
            ],
      ext_modules= [
            Extension("pygmin.mindist.overlap", 
                      ["pygmin/mindist/overlap.f90"]),
                    
            Extension("pygmin.optimize.mylbfgs_fort", 
                      ["pygmin/optimize/mylbfgs_fort.f90"]),
            Extension("pygmin.optimize.mylbfgs_updatestep", 
                      ["pygmin/optimize/mylbfgs_updatestep.f90"]),
                    
            Extension("pygmin.potentials.fortran.AT", 
                      ["pygmin/potentials/fortran/AT.f"]),
            #Extension("pygmin.potentials.fortran.lj_hess", 
            #          ["pygmin/potentials/fortran/lj_hess.f"]),
            Extension("pygmin.potentials.fortran.ljpshiftfort", 
                      ["pygmin/potentials/fortran/ljpshiftfort.f"]),
            Extension("pygmin.potentials.fortran.lj", 
                      ["pygmin/potentials/fortran/lj.f90"]),
            Extension("pygmin.potentials.fortran.ljcut", 
                      ["pygmin/potentials/fortran/ljcut.f90"]),
            Extension("pygmin.potentials.fortran.rmdrvt", 
                      ["pygmin/potentials/fortran/rmdrvt.f90"]),
            Extension("pygmin.potentials.fortran.soft_sphere_pot", 
                      ["pygmin/potentials/fortran/soft_sphere_pot.f90"]),
                    
            Extension("pygmin.utils._fortran_utils", 
                      ["pygmin/utils/_fortran_utils.f90"]),
                    
            Extension("pygmin.potentials.cpp.lj", 
                      ["pygmin/potentials/cpp/lj.cpp"],
                      include_dirs=["include"]),
            Extension("pygmin.potentials.stock.stockaa", 
                      ["pygmin/potentials/stock/stockaa.cpp"],
                      include_dirs=["include"]),
            Extension("pygmin.potentials.salt_.matrix", 
                      ["pygmin/potentials/salt_/matrix.cpp"],
                      include_dirs=["include"]),
            Extension("pygmin.potentials.salt_.salt", 
                      ["pygmin/potentials/salt_/salt.cpp"],
                      include_dirs=["include"]),
                    
                ]
        )

