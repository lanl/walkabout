module disptensor_class

use precision 
use isotropictensor_class 
use lichtnertensor_class 
use bftensor_class 

implicit none 

private 

public disptensor 
public disptensor_ 
public calculateBtensor 
public calculateDtensor 
public MaxTransverseDispersion 

type disptensor 
  type(isotropictensor) :: isotropic 
  type(lichtnertensor) :: lichtner 
  type(bftensor) :: burnettfrind 
  integer :: type  
end type disptensor 

integer, parameter :: LICHT=1
integer, parameter :: ISOTR=0
integer, parameter :: BF=2

real(kind=dp) :: dm = 8.64d-5 ! molecular diffusion m2/day

contains 

function disptensor_( funit) result(aobj) 

type(disptensor) :: aobj 
character(len=5) :: option 
integer :: funit 


read(funit,'(a)') option 

select case( option ) 
case('ISOTR') 
 aobj%isotropic = isotropictensor_(funit) 
 aobj%type=ISOTR
case('LICHT') 
 aobj%lichtner = lichtnertensor_(funit) 
 aobj%type=LICHT
case('BF') 
 aobj%burnettfrind = bftensor_(funit) 
 aobj%type=BF
case DEFAULT 
  print *, 'tensor option not defined ', option 
end select 

end function disptensor_ 

!*********************************************************************

function calculateBtensor(aobj , vel ) result(amat) 

type(disptensor) :: aobj 
real(kind=dp), dimension(3,3) :: amat 
real(kind=dp), dimension(3) :: vel 
 
if( aobj%type .eq. ISOTR) then 
  amat=cB( aobj%isotropic, vel ) 
else if( aobj%type .eq. LICHT) then 
  amat=cB( aobj%lichtner, vel ) 
else if( aobj%type .eq. BF) then 
  amat=cB( aobj%burnettfrind, vel ) 
else 
  print *, ' tensor object not defined'  
  stop 
end if 

return

end function calculateBtensor
!*********************************************************************
function calculateDtensor(aobj , vel ) result(amat) 

type(disptensor) :: aobj 
real(kind=dp), dimension(3,3) :: amat 
real(kind=dp), dimension(3) :: vel 
 
if( aobj%type .eq. ISOTR) then 
  amat=cD( aobj%isotropic, vel ) 
else if( aobj%type .eq. LICHT) then 
  amat=cD( aobj%lichtner, vel ) 
else if( aobj%type .eq. BF) then 
  amat=cD( aobj%burnettfrind, vel ) 
else 
  print *, ' tensor object not defined'  
  stop 
end if 

return

end function calculateDtensor

!*********************************************************************

subroutine MaxTransverseDispersion( aobj , icell, alphaT, Dm) 

type(disptensor) :: aobj 
integer :: icell 
real(kind=dp) :: alphaT, Dm 
 
if( aobj%type .eq. ISOTR) then 
  alphaT= mTD( aobj%isotropic, icell) 
else if( aobj%type .eq. LICHT) then 
  alphaT= mTD( aobj%lichtner, icell) 
else if( aobj%type .eq. BF) then 
  alphaT= mTD( aobj%burnettfrind, icell) 
  Dm = aobj%burnettfrind%Dm  
else 
  print *, ' tensor object not defined'  
  stop 
end if 

end subroutine MaxTransverseDispersion 


end module disptensor_class 
