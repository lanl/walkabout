module linked_list 

implicit none 

private 

public add2queue
public add2queuend 
public fromqueue  
public destroyqueue  

type record 
 integer :: i
 type(record), pointer :: next => null() 
end type record 

type(record), pointer :: list,  first , current 

integer :: count 


contains 

!********************************************

subroutine add2queue(i) 
integer :: i 
logical :: empty 


empty=.true. 
if(  associated(first) ) empty = .false. 
if( empty ) then  
  allocate(list)  
  nullify(list%next) 
  current => list
  first => list 
  count = 0 
end if  

list%i=i 

allocate(list%next)  
list=>list%next; 
nullify(list%next) 

count = count + 1

if( empty ) current%next => list  

end subroutine add2queue 

!*******************************************

function fromqueue() result(i)  

implicit none 

integer :: i,istat 
logical :: test 
type(record),pointer :: previous 


 if( .not. associated(current%next) ) then
   i=-1 
   return 
 end if  

 i=current%i 
 previous =>  current 
 current => current%next 

! deallocate(previous, stat=istat ) 
! if( istat .ne. 0) then 
!   print *, istat 
!   print *, i 
!   stop 
! end if 

! count = count -1 

 return  

end function fromqueue  

!*******************************************


subroutine add2queuend(i) 
integer :: i , j 
logical :: empty 
type(record),pointer :: ptr 


empty=.true. 
if(  associated(first) ) empty = .false. 
if( empty ) then  
  allocate(list)  
  nullify(list%next) 
  current => list
  first => list 
  count = 0 
end if  

!ptr => first  
!do j=1,count 
!   if( ptr%i == i) return 
!   ptr => ptr%next
!end do 

! if here then element was not already in list 

list%i=i 

allocate(list%next)  
list=>list%next; 
nullify(list%next) 

count = count + 1

if( empty ) then 
   current%next => list  
   first%next => list 
end if 

end subroutine add2queuend   

!*******************************************

function destroyqueue() result(i)  

implicit none 

integer :: i,istat 
logical :: test 
type(record),pointer :: previous 


 if( count == 0  ) then
   i=-1 
   deallocate(first) 
   return 
 end if  


 i=first%i 
 previous =>  first  
 first => first%next 


 deallocate(previous, stat=istat ) 
 if( istat .ne. 0) then 
   print *, istat 
   print *, i 
   stop 
 end if 

 count = count -1 

 return  

end function destroyqueue

end module linked_list  
