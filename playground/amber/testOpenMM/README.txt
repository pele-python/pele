OpenMM 6 installation from source on Ubuntu 13.10
=================================================
Download OpenMM6.0-Linux64.zip fromi Downloads page on https://simtk.org/home/openmm

$ cd OpenMM6.0-Source 
$ mkdir build_openmm
$ cd  build_openmm
$ ccmake -i .. 
# change CMAKE_INSTALL_PREFIX to ~/opt/bin/openmm6 
$ make install 
$ sudo make PythonInstall 

Installed files: 
$ ls ~/opt/bin/openmm6/
bin  docs  examples  include  lib  licenses
$ ls /usr/local/lib/python2.7/dist-packages/
OpenMM-5.2.0.egg-info  OpenMM-6.0.0.egg-info  simtk

Test
=====

$ cd ~/opt/bin/openmm6/examples
~/opt/bin/openmm6/examples$ python testInstallation.py 
There are 1 Platforms available:

1 Reference - Successfully computed forces

