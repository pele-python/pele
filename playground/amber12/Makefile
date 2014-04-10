mypotential.so: mypotential.o _mypotential.o
	gfortran -Wall -Wall -shared mypotential.o  _mypotential.o -lgfortran -o mypotential.so

mytest: mytest.o _mypotential.o
	gcc -Wall mytest.o _mypotential.o -o mytest -lgfortran

mytest.o: mytest.c
	gcc -Wall -g -std=c99 -c mytest.c

mypotential.o: mypotential.c
	gcc -pthread -fno-strict-aliasing -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -c mypotential.c  `python-config --include` -o mypotential.o

mypotential.c: mypotential.pyx
	cython mypotential.pyx


_mypotential.o: _mypotential.f90
	gfortran -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -c _mypotential.f90 -o _mypotential.o
