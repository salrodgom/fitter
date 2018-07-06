install:
	gfortran fitter.f90 -o fitter -ffree-form -march=native
all:
	gfortran fitter.f90 -o fitter -ffree-form -march=native
	cp lib/peros.lib.original peros.lib ; ./fitter < input
clean:
	rm -rf fitter *.o *.mod
