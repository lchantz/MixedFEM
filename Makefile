all: compile run

fortran_files := $(shell ls | grep "\.f90")
compile: 
	gfortran $(fortran_files) \
	-o MixedFEM \
	-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none \
	-fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow \
	-finit-real=nan \
	# ./lapack-3.11/libtmglib.a \
	# ./lapack-3.11/liblapacke.a \
	# ./lapack-3.11/liblapack.a \
	# ./lapack-3.11/librefblas.a \

run: 
	./MixedFEM