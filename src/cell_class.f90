module cell_class 

use precision  
implicit none 

private 
public cell 

type cell 
 real(kind=dp), dimension(3) :: posi 
 integer :: nn 
 integer, dimension(:), pointer :: nbrs 
 integer, dimension(:), pointer :: indx2area 
 integer, dimension(:), pointer :: indx2conn 
 integer :: NumNoflow 
 logical :: inflow 
 real(kind=dp), pointer, dimension(:,:) :: bmat !to be dimensioned NumNoflow,3 
end type cell 

end module cell_class 

