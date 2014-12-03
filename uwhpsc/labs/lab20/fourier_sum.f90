
program fourier_sum

    implicit none
    integer :: n, i, nterms, k
    real(kind=8) :: pi, ck, term, dx
    real(kind=8), dimension(:), allocatable :: x, fsum
    open(unit=21, file='input_data.txt', status='old')

    read(21,*) n
    read(21,*) nterms

    print '("n = ",i6)', n
    print '("nterms = ",i6)', nterms

    allocate(x(0:n+1), fsum(0:n+1))
    pi = acos(-1.d0)

    open(unit=22, file='xf.txt', status='unknown')
    open(unit=23, file='frames.txt', status='unknown')

    dx = pi / (n + 1.d0)

    do i=0,n+1
        x(i) = i*dx
        write(22,221) x(i), exp(x(i))
221     format(2e20.10)
        fsum(i) = 0.d0
        enddo

    do k=1,nterms
        ck = 2*k*(1. - (-1.)**k * exp(pi)) / ((k**2 + 1)*pi)
        do i=0,n+1
            term = ck * sin(k*x(i))
            fsum(i) = fsum(i) + term
            write(23,221) term, fsum(i)
            enddo
        enddo

end program fourier_sum
