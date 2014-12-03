
Errors in array_omp.f90:


Missing 
    use omp_lib

Calling omp_get_num_threads outside a parallel block will always return 1
since only the master thread is active.

The statment
    $omp parallel do private(i)
has two errors:
    should be !$omp not $omp
    private(i) should not be declared for parallel do j loop inside do i loop.

