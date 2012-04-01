all: fortran

fortran:
	cd potentials/fortran ; make

cpp:
	cd potentials/cpp; cmake . ; make

