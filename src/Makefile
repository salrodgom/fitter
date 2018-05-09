# -fopenmp
install:
	gfortran fitter.f90 -o fitter
all:
	gfortran fitter.f90 -o fitter
	cp lib/peros.lib.original peros.lib ; ./fitter < input
clean:
	rm -rf fitter *.o *.mod
