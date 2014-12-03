
FC = gfortran
FFLAGS = -O3
LFLAGS = -lblas -llapack
.PHONY: clean test plot

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

hilbert.exe: hilbert.o 
	$(FC) hilbert.o  $(LFLAGS) -o hilbert.exe

hilbert2.exe: hilbert2.o 
	$(FC) hilbert2.o  $(LFLAGS) -o hilbert2.exe

test: hilbert.exe
	./hilbert.exe

cond.txt: hilbert2.exe
	./hilbert2.exe

plot: cond.txt
	python plot_cond.py

clean:
	rm -f *.o *.exe
