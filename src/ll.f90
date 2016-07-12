module linked_list 

implicit none 

private 

!public add2queue
public add2queuend 
public fromqueue  
public destroyqueue  

type record 
 integer :: i
 type(record), pointer :: next => null() 
end type record 

type, public :: queue
 integer :: count
 type(record), pointer :: list,current,first => null()
end type queue

contains 

!********************************************

!subroutine add2queue(i) 
!integer :: i 
!logical :: empty 


!empty=.true. 
!if(  associated(first) ) empty = .false. 
!if( empty ) then  
 ! allocate(list)  
  !nullify(list%next) 
  !current => list
  !first => list 
  !count = 0 
!end if  

!list%i=i 

!allocate(list%next)  
!list=>list%next; 
!nullify(list%next) 

!count = count + 1

!if( empty ) current%next => list  

!end subroutine add2queue 

!*******************************************

function fromqueue(q) result(i)  
implicit none 

type(queue) :: q
integer :: i,istat 
logical :: test 
type(record),pointer :: previous 


 if( .not. associated(q%current%next) ) then
   i=-1 
   return 
 end if  

 !q%current = q%list
 i=q%current%i 
 previous =>  q%current 
 q%current => q%current%next 

 !deallocate(previous, stat=istat ) 
 !if( istat .ne. 0) then 
  ! print *, istat 
   !print *, i 
   !stop 
 !end if 

!q%count = q%count -1 

 return  

end function fromqueue  

!*******************************************


subroutine add2queuend(q,i)
type(queue) :: q 
integer :: i , j 
logical :: empty 
type(record),pointer :: ptr 


empty=.true. 
if(  associated(q%first) ) empty = .false. 
!nullify(q%list)
if( empty ) then  
  allocate(q%list)  
  nullify(q%list%next) 
  !q%list%i=i
  q%current => q%list
  q%first => q%list 
  q%count = 0 
end if  

!q%current => q%list
!q%first => q%list 
!q%count = 0 

ptr => q%first  
do j=1,q%count 
   if( ptr%i == i) return 
   ptr => ptr%next
end do 

! if here then element was not already in list 

q%list%i=i 

allocate(q%list%next)  
!allocate(q%current%next)
q%list=>q%list%next; 
nullify(q%list%next) 


q%count = q%count + 1

if( empty ) then 
   !q%current => q%list
   !q%first => q%list
   !q%count = 0
   q%current%next => q%list  
   q%first%next => q%list 
end if 

!q%count = q%count + 1

!if(  associated(q%first) ) empty = .false.

end subroutine add2queuend   

!*******************************************

function destroyqueue(q) result(i)  

implicit none 

type(queue) :: q
integer :: i,istat 
logical :: test 
type(record),pointer :: previous 


 if(q%count == 0  ) then
   i=-1 
   deallocate(q%first) 
   return 
 end if  


 i=q%first%i 
 previous =>  q%first  
 q%first => q%first%next 


 deallocate(previous, stat=istat ) 
 if( istat .ne. 0) then 
   print *, istat 
   print *, i 
   stop 
 end if 

 q%count = q%count -1 

 return  

end function destroyqueue

end module linked_list  
