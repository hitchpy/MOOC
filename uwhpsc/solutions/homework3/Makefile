
OBJECTS = random_util.o approx_norm.o main.o
FFLAGS = -O2 -fopenmp
#FFLAGS = -g -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid 

.PHONY: data test clean 

# default values for input_data.txt:

n = 50
seed = 12345
nsamples = 1000
nthreads = 2
method = 1

data:
	@echo $(n) > input_data.txt
	@echo $(seed) >> input_data.txt
	@echo $(nsamples) >> input_data.txt
	@echo $(nthreads) >> input_data.txt
	@echo $(method) >> input_data.txt
	@echo Wrote data to input_data.txt


test: test.exe data
	./test.exe

test.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

