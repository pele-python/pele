from distutils.core import setup

setup(name='PyGMIN', 
        version='0.1', 
        description="Pythin implementation of GMIN, OPTIM, and PATHSAMPLE",
        url='http://github.com/js850/minGMIN.git',
        packages=["pygmin", 
            "pygmin.potentials",
            "pygmin.gui",
            "pygmin.mindist",
            "pygmin.optimize",
            "pygmin.printing",
            "pygmin.takestep",
            "pygmin.utils",
            "pygmin.wham"
            ]
        )

