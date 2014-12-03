

To compile the version with parallel do loops:
    $ gfortran -fopenmp random_util.f90 mean.f90

To compile the coarse grain parallel version:
    $ gfortran -fopenmp random_util.f90 mean2.f90
