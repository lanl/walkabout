module bftensor_class

use precision 

implicit none 

private 

public bftensor 
public bftensor_ 
public cD 
public cB 
public mTD 

type bftensor 
 real(kind=dp) :: alphaL , alphaTv, alphaTh, Dm 
end type bftensor 

interface cB 
  module procedure cB_bf 
end interface cB 

interface cD 
  module procedure cD_bf 
end interface cD 

interface mTD  
  module procedure mTD_bf 
end interface mTD 


contains 

!*********************************************************************************

function bftensor_(funit ) result(aobj)  
  type(bftensor), target :: aobj 
  integer :: funit 

  read(funit,*) aobj%alphaL, aobj%alphaTh, aobj%alphaTv , aobj%Dm 


end function bftensor_ 

!********************************************************************************
function cB_bf(aobj,vel) result(disp) 

  type(bftensor) :: aobj 
  real(kind=dp), dimension(3) :: vel 
  real(kind=dp), dimension(3,3) :: disp  
  real(kind=dp) :: alphaL, alphaTh, alphaTv, absvel 
  real(kind=dp) :: group1, group2, group3  , foo2 
  real(kind=dp) :: Dm 

!  real(kind=dp), parameter :: Dm = 8.64d-5 ! molecular diffusion m2/day

  alphaL=aobj%alphaL 
  alphaTh=aobj%alphaTh 
  alphaTv=aobj%alphaTv 
  Dm=aobj%Dm 

  ! displacement tensor (not dispersion tensor)  

  if( alphaL + alphaTh + alphaTv + Dm .le. 1.0d-20) then 
   disp = 0.0d0 
   return 
  end if 

  absvel=sqrt(dot_product(vel,vel)) + 1.0d-20 
  foo2=vel(1)*vel(1)+vel(2)*vel(2) + 1.0d-20 
  group1 = sqrt(2.0d0*(alphaL*absvel + Dm) )
  group2 = sqrt(2.0d0*(alphaTv*absvel + Dm) )/sqrt(foo2) 
  group3 = sqrt(2.0d0* ( ( alphaTh*foo2+alphaTv*vel(3)*vel(3) )/absvel + Dm )) 

  disp(1,1)= vel(1)*group1 / absvel  
  disp(1,2)= -vel(1)*vel(3)*group2 / absvel  
  disp(1,3)= -vel(2)/sqrt(foo2) * group3  
  disp(2,1)= vel(2)*group1 / absvel  
  disp(2,2)= -vel(2)*vel(3)*group2 / absvel  
  disp(2,3)= vel(1)/sqrt(foo2) * group3  
  disp(3,1)= vel(3)*group1 / absvel  
  disp(3,2)= sqrt(2.0*foo2/(absvel*absvel) * ( alphaTv*absvel + Dm) ) 
  disp(3,3)=0.0d0 



end function cB_bf

!********************************************************************************
function cD_bf(aobj,vel) result(dcoef) 

  type(bftensor) :: aobj 
  real(kind=dp), dimension(3) :: vel 
  real(kind=dp), dimension(3,3) :: dcoef 
  real(kind=dp) :: alphaL, alphaTh, alphaTv, absvel 
  real(kind=dp) :: group1, group2, group3  , foo2 

  alphaL=aobj%alphaL 
  alphaTh=aobj%alphaTh 
  alphaTv=aobj%alphaTv 

  
  ! dispersion tensor  

  absvel=sqrt(dot_product(vel,vel)) + 1.0d-20 
  foo2=vel(1)*vel(1)+vel(2)*vel(2) 

  dcoef(1,1)= alphaL*vel(1)*vel(1) + alphaTh*vel(2)*vel(2) + alphaTv*vel(3)*vel(3)  
  dcoef(1,2)= (alphaL-alphaTh)*vel(2)*vel(1)  
  dcoef(1,3)= (alphaL-alphaTv)*vel(1)*vel(3)   
  dcoef(2,1)= dcoef(1,2)  
  dcoef(2,2)= alphaTh*vel(1)*vel(1) + alphaL*vel(2)*vel(2) + alphaTv*vel(3)*vel(3) 
  dcoef(2,3)= (alphaL-alphaTv)*vel(2)*vel(3)  
  dcoef(3,1)= dcoef(1,3)  
  dcoef(3,2)= dcoef(2,3)  
  dcoef(3,3)= alphaTv*foo2 + alphaL*vel(3)*vel(3) 

  dcoef=dcoef/absvel 


end function cD_bf

!*******************************************************************************
function mTD_bf( aobj, icell) result(aval) 

 type(bftensor) :: aobj
 integer :: icell 
 real(kind=dp ) :: aval

aval = max( aobj%alphaTv, aobj%alphaTh, aobj%alphaL) 

end function  mTD_bf


end module bftensor_class 

