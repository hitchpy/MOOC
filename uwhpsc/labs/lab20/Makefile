
OBJECTS = fourier_sum.f90

.PHONY: data test clean movie clobber

# default values for input_data.txt:

n = 200
nterms = 20


data:
	@echo $(n) > input_data.txt
	@echo $(nterms) >> input_data.txt

test: test.exe data
	./test.exe

test.exe: $(OBJECTS)
	gfortran $(OBJECTS) $(LFLAGS) -o test.exe
    
%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

movie: test
	python animate.py

clean:
	rm -f *.o *.exe *.mod *.pyc

clobber: clean
	rm -f xf.txt frames.txt FourierSum.html
	rm -rf _plots
