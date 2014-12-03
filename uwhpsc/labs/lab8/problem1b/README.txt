
The errors in the original:

In the Makefile, constants.o must come before main.o since main uses this
module.  If constants.mod already exists it will still compile, but
following "make clean" it will not compile as written.

constants.f90 is missing a "contains" statement before the subroutine.

area.f90: the variable a must be intent(out)


