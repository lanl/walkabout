module isotropictensor_class 

use precision 

implicit none 

private 

public isotropictensor 
public isotropictensor_ 
public cB
public cD
public mTD 

type isotropictensor 
 real(kind=dp) :: alphaL , alphaT 
end type isotropictensor 

interface cB 
  module procedure cB_iso
end interface cB 

interface cD 
  module procedure cD_iso
end interface cD 

interface mTD 
  module procedure mTD_iso
end interface mTD 

contains 

!*********************************************************************************

function isotropictensor_(funit ) result(aobj)  
  type(isotropictensor), target :: aobj 
  integer :: funit 

  read(funit,*) aobj%alphaL, aobj%alphaT


end function isotropictensor_ 

!********************************************************************************
function cD_iso(aobj,vel) result(dcoef) 

  type(isotropictensor) :: aobj 
  real(kind=dp), dimension(3) :: vel 
  real(kind=dp), dimension(3,3) :: dcoef  
  real(kind=dp) :: alphaL, alphaT , absvel 

  alphaL=aobj%alphaL 
  alphaT=aobj%alphaT 

  
  ! dispersion coefficient
  absvel=sqrt(dot_product(vel,vel))
  dcoef(1,1)=alphaT*absvel + (alphaL-alphaT)* vel(1)*vel(1)/absvel
  dcoef(1,2)= (alphaL-alphaT)* vel(1)*vel(2)/absvel
  dcoef(1,3)= (alphaL-alphaT)* vel(1)*vel(3)/absvel
  dcoef(2,1)= (alphaL-alphaT)* vel(2)*vel(1)/absvel
  dcoef(2,2)=alphaT*absvel + (alphaL-alphaT)* vel(2)*vel(2)/absvel
  dcoef(2,3)= (alphaL-alphaT)* vel(2)*vel(3)/absvel
  dcoef(3,1)= (alphaL-alphaT)* vel(3)*vel(1)/absvel
  dcoef(3,2)= (alphaL-alphaT)* vel(3)*vel(2)/absvel
  dcoef(3,3)=alphaT*absvel + (alphaL-alphaT)* vel(3)*vel(3)/absvel

end function cD_iso

!********************************************************************************
function cB_iso(aobj,vel) result(disp) 

  type(isotropictensor) :: aobj 
  real(kind=dp), dimension(3) :: vel 
  real(kind=dp), dimension(3,3) :: disp
  real(kind=dp) :: alphaL, alphaT , absvel 

  alphaL=aobj%alphaL 
  alphaT=aobj%alphaT 

  ! dispersion coefficient
  absvel=sqrt(dot_product(vel,vel))

   print *, 'displacement not implemented in isotropic tensor' 
   stop 

end function cB_iso

!*******************************************************************************
function mTD_iso( aobj, icell) result(aval) 

 type(isotropictensor) :: aobj 
 integer :: icell 
 real(kind=dp ) :: aval

 aval = aobj%alphaT 

end function mTD_iso 

end module isotropictensor_class 

