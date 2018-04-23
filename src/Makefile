# -fopenmp
all:
	gfortran fitter.f90 -o fitter -fopenmp
	cp lib/peros.lib.original peros.lib ; ./fitter < input
clean:
	rm -rf fitter *.o *.mod
