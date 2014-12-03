
program demo

! To demonstrate trapping on divide by zero and use of gdb.
! Similar to lab7/demo.py.


real(kind=8) :: x(1001),y(1001)
integer :: j

do j=1,1001
    x(j) = -20.d0 + (j-1)*40.d0/1000.d0
    y(j) = sin(x(j)) / x(j)
    enddo
        
print *, "The max value of y is ",maxval(y)

print *, "x(501) is ", x(501)
print *, "y(501) is ", y(501)
    
end program demo
