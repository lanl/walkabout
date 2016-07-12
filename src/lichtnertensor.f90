module lichtnertensor_class 

use precision 

implicit none 

private 

public lichtnertensor 
public lichtnertensor_ 
public cB
public cD
public mTD 

type lichtnertensor 
 real(kind=dp) :: alphaL , alphaT 
end type lichtnertensor 

interface cB 
  module procedure cB_licht
end interface cB 

interface cD 
  module procedure cD_licht
end interface cD 

interface mTD 
  module procedure mTD_licht
end interface mTD 

contains 

!*********************************************************************************

function lichtnertensor_(funit) result(aobj)  
  type(lichtnertensor)  :: aobj 
  integer :: funit 

  read(funit,*) aobj%alphaL, aobj%alphaT

end function lichtnertensor_ 

!********************************************************************************
function cD_licht(aobj,vel) result(dcoef) 

  type(lichtnertensor) :: aobj 
  real(kind=dp), dimension(3) :: vel 
  real(kind=dp), dimension(3,3) :: dcoef  
  real(kind=dp) :: alphaL, alphaT , absvel 
  
  alphaT=aobj%alphaT 
  alphaL=aobj%alphaL 

  ! dispersion coefficient
  absvel=sqrt(dot_product(vel,vel))

  print *, 'displacement not implemented in licthtner tensor' 
  stop 

end function cD_licht 

function cB_licht(aobj,vel) result(disp) 

  type(lichtnertensor) :: aobj
  real(kind=dp), dimension(3) :: vel
  real(kind=dp), dimension(3,3) :: disp
  real(kind=dp) :: alphaL, alphaT , absvel

  alphaL=aobj%alphaL
  alphaT=aobj%alphaT

  ! dispersion coefficient
  absvel=sqrt(dot_product(vel,vel))

  print *, 'displacement not implemented in lichtner tensor' 
  stop 

end function cB_licht 

!*******************************************************************************
function mTD_licht( aobj, icell) result(aval) 

 type(lichtnertensor) :: aobj
 integer :: icell
 real(kind=dp ) :: aval

 aval = aobj%alphaT 
! aval = max( aobj%alphaTv, aobj%alphaTl) 

end function mTD_licht


end module lichtnertensor_class 

