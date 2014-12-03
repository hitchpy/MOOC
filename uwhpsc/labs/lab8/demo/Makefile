
DEBUG_FLAGS = -g -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid 
.PHONY: run debug

run:
	gfortran demo.f90
	./a.out

debug:
	gfortran $(DEBUG_FLAGS) demo.f90
	./a.out


