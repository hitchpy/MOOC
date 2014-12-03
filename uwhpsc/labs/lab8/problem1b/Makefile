

OBJECTS = constants.o main.o area.o 
.PHONY: run, clean

run: main.exe
	./main.exe

output.txt: main.exe
	./main.exe > output.txt

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe

%.o : %.f90
	gfortran -c $< 

clean:
	rm -f $(OBJECTS) main.exe constants.mod

