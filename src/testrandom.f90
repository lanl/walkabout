program testrandom 

implicit none 

integer :: i 
integer, dimension(8) :: seed 
integer, parameter :: dp=kind(1.0d0) 
real(kind=dp), parameter :: ascale=sqrt(3.0d0) 
real(kind=dp) :: r 

call system_clock(seed(1))
call random_seed(put=seed)

do i=1,1000 
call random_number(r) 
r=(2.0d0*r-1.0d0) * ascale 
print *, r 

end do 


end program testrandom 

