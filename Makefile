FC=gfortran
FFLAGS= -ffixed-form -ffixed-line-length-none -g -Wall

utils=ljpshiftfort

all: ${utils}.so

${utils}.o: ${utils}.f90
	$(FC) $(FFLAGS) -c -o $@ $<

${utils}.so: ${utils}.f
	f2py -c -m ${utils} --fcompiler=gfortran  --link-lapack_opt $< > ${utils}.setup

#f2py --opt="-O3" -c -m fd_rrt1d --fcompiler=gfortran  --link-lapack_opt *.f
