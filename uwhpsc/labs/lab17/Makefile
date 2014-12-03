
FC = gfortran
FFLAGS = -O3
LFLAGS = -lblas -llapack
.PHONY: clean test plot

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

tridiag.exe: tridiag.o 
	$(FC) tridiag.o  $(LFLAGS) -o tridiag.exe

test: tridiag.exe
	./tridiag.exe

clean:
	rm -f *.o *.exe
