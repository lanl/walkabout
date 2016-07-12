module iotools 

implicit none 

public findunitnum 

contains 

!**********************************************************************
 function findunitnum() result(iunit)

  integer :: iunit
  logical :: op
  integer,parameter :: maxunit=99 
  

  do iunit=10,maxunit
    inquire(unit=iunit, opened = op)
    if( .not. op )return
  end do

  print *, ' all file unit numbers used'
  stop

 end function findunitnum 

subroutine printheader( funit, ptitle ) 
integer :: funit 
character(len=8) :: date
character(len=11) :: time 
character(len=30) :: vernum="     walkabout 1.2             " 
character(len=5) :: zone 
character(len=80) ptitle 

call date_and_time( date, time, zone) 

write(funit,*) vernum, date, time 
write(funit,*) ptitle 

end subroutine printheader 

subroutine printversion 

character(len=40) name 

name =   'Walkabout 1.2 November 2011'  

print *, '                              ' 
print *, name 
print *, '        Scott Painter         ' 
print *, 'Los Alamos National Laboratory'  
print *, '                              ' 

end subroutine printversion 

end module iotools  
