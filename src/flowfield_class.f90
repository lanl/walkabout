module flowfield_class

use precision 
use iotools 
use mesh_class 
use cell_class 
use tracktools , only : invert3x3, invert2x2 
use disptensor_class 

implicit none 

private 
public flowfield, flowfield_, reconstruct_flowfield 


type flowfield  
 real(kind=dp) :: time  ! days
 real(kind=dp), dimension(:), pointer :: fluxes ! initially in kg/s then converted to m3/day 
 real(kind=dp), dimension(:,:), pointer :: vfield ! m/day 
 real(kind=dp), dimension(:,:,:), pointer :: disp ! displacement matrix 
 real(kind=dp), dimension(:,:,:), pointer :: dcoef ! dcoef matrix 
 real(kind=dp), dimension(:), pointer :: dtdisp ! time step limit by dispersive 
 real(kind=dp), dimension(:), pointer :: porosity 
 real(kind=dp), dimension(:), pointer :: saturation  
 real(kind=dp), dimension(:), pointer :: density ! kg/m3  
 type(disptensor), dimension(:), pointer  :: tensor ! can make pointer array to save space
 type(mesh), pointer :: ptr2mesh   ! allows for adaptive mesh 
 logical amaread
end type flowfield  

contains 

!*********************************************************************
function flowfield_(namefile,  ptr2mesh) result( aobj) 

! loads the flowfield object 

character(len=*) :: namefile  
character(len=70) :: controlfile, finfile , avsfile , zonefile , amafile
character(len=80) :: line  
character(len=80) :: buffer  
type(flowfield) :: aobj 
type(mesh), pointer :: ptr2mesh 
type(disptensor) :: atensor  
integer :: nvals 
integer :: funit1, funit2 
integer :: isat,ipor, iden, ncol 
integer :: i, i1, i2, i3 
real(kind=dp), dimension(:), allocatable :: tmparr 
integer, dimension(:), allocatable :: itmparr 
integer :: ios 
integer :: azone, bzone , nnum
logical finread, amaread

finread = .false.
amaread = .false.

aobj%amaread = .false.


! liquid fluxes are in the FEHM fin file 
! initially in kg/s 
funit1=findunitnum() 
open(unit=funit1,file=namefile,action='read')
 do
  read(unit=funit1,fmt='(a)', end=994) buffer
  buffer=adjustl(buffer)
  if( buffer(1:4) .eq. "fin:") exit
 end do
 finfile = buffer(5:)
 finread = .true.
 994 continue
 close(funit1)
 
funit1=findunitnum() 
open(unit=funit1,file=namefile,action='read')
 do
  read(unit=funit1,fmt='(a)', end=993) buffer
  buffer=adjustl(buffer)
  if( buffer(1:4) .eq. "ama:") exit
 end do
 amafile = buffer(5:)
 amaread = .true.
 aobj%amaread = .true.
 993 continue
 close(funit1)
 

if(finread .and. amaread)then
	print *, "Both FEHM and Amanzi input files specified.  Please specify one or the other."
	stop
else if(.not. finread .and. .not. amaread)then
	print *, "Neither FEHM or Amanzi input files specified.  Please specify one or the other."
	stop
endif

if(finread)then
	funit1=findunitnum() 
	open(unit=funit1,file=finfile,action='read')
	do 
 		read(funit1,'(a80)') line  
 		if( line .eq. 'liquid flux') exit 
	end do 

	read(funit1,*) nvals 
	allocate(aobj%fluxes(nvals) ) 
	read(funit1,*) aobj%fluxes  

	! convert to kg/day 
	aobj%fluxes = aobj%fluxes * 3600.0d0*24.0d0

	close(funit1) 
endif

if(amaread)then
	funit1=findunitnum() 
	open(unit=funit1,file=amafile,action='read')
	read(funit1, *) nvals
	close(funit1)
endif

aobj%ptr2mesh => ptr2mesh 

!nvals = ptr2mesh%ncells 

print *, 'nvals are', nvals
allocate( aobj%vfield(nvals,3) ) 
allocate( aobj%disp(nvals,3,3) )  
allocate( aobj%dcoef(nvals,3,3) )  
allocate( aobj%dtdisp(nvals) ) 
allocate( aobj%porosity(nvals) ) 
allocate( aobj%density(nvals) ) 
allocate( aobj%saturation(nvals) ) 
allocate( aobj%tensor(nvals) ) 

if(amaread)then
	funit1=findunitnum() 
	open(unit=funit1,file=amafile,action='read')
	!print *, nvals
        read(funit1,*)  ! skip the first line for nvals
	do i=1,nvals
		read(funit1, *), buffer, aobj%vfield(i,1), aobj%vfield(i,2), aobj%vfield(i,3)
	enddo
	close(funit1)
        aobj%vfield = aobj%vfield*3600.0*24.0   ! converted to m/day
endif

aobj%time = 0.0d0 

 ! search the namefile to get name of control file 
 funit1=findunitnum() 
 open(unit=funit1,file=namefile,action='read')

 do
  read(unit=funit1,fmt='(a)', end=999) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8) .eq. "control:") exit
 end do
 controlfile = buffer(9:)
 close(funit1)

! *************************************************
! dispersion tensor block 
 funit1=findunitnum() 
 open(unit=funit1, file=controlfile, action='read') 
 do
  read(unit=funit1,fmt='(a)', end=998) buffer
  buffer=adjustl(buffer)
  if( buffer(1:7).eq.'DTENSOR') exit  
 end do

 ! loop over all zone defs 
 do 

  read(funit1,'(a)') buffer
  buffer=adjustl(buffer)
  if(buffer(1:3) .eq. 'END') exit  

  read(unit=buffer, fmt=*, iostat=ios) i1, i2, i3

  atensor=disptensor_(funit1) ! read then assign below 

  if( ios .ne. 0 ) then 

   ! line contains a file name and a zone name 

    read(buffer,*) zonefile, azone 
    funit2=findunitnum() 
    open(unit=funit2, file=zonefile, action='read') 
    
    read(funit2,'(a)') buffer
    buffer=adjustl(buffer) 
    if( buffer(1:4) .ne. 'zone') goto 995  

    do 
992  read(unit=funit2, fmt=*, end=995) bzone  
     if( bzone .eq. azone) then
        exit 
     else
        read(unit=funit2,fmt=*, end=995) buffer
        read(unit=funit2,fmt=*, end=995) nnum
        allocate( itmparr(nnum) ) 
        read(unit=funit2, fmt=*, end=995) itmparr 
        deallocate(itmparr)
        goto 992     
     endif
    end do 

    ! if here then zone was found 
    read(funit2,'(a)') buffer
    buffer=adjustl(buffer) 
    if(buffer(1:4) .ne. 'nnum') goto 995 

    read(unit=funit2, fmt=*, end=995) nnum  
    allocate( itmparr(nnum) ) 
    read(unit=funit2, fmt=*, end=995) itmparr 
    close(funit2)      ! added by Zhiming. Otherwise cannot open for the next zone, if all zones are in one file.
    do i=1,nnum 
     aobj%tensor( itmparr(i) ) = atensor  
    end do 
    deallocate(itmparr) 

  else  ! assign by cell numbers 

    if(i1.eq. 1 .and. i2.eq. 0 .and. i3.eq. 0) then 
     i1=1 
     i2=nvals 
     i3=1 
    end if 

    do i=i1, i2, i3 
     aobj%tensor(i) = atensor  
    end do 

  end if 
 
 end do 

 ! now check that all cells have been assigned a disptensor 
 do i=1,nvals 
   if( aobj%tensor(i)%type .lt. 0) then 
     print *, 'cell ', i, ' has not been assigned a dispersion tensor. stopping' 
     stop 
   end if 
 end do 

close(funit1) 
! end of dispersion tensor block 
! *******************************************************************************

 ! search the namefile to get name of avs file
 ! used for porosity and saturation 
 
 if( .not. amaread)then
 	funit1=findunitnum()
	 open(unit=funit1,file=namefile,action='read')
	 do
	  read(unit=funit1,fmt='(a)', end=997) buffer
	  buffer=adjustl(buffer)
	  if( buffer(1:4) .eq. "avs:") exit
	 end do
	 avsfile = buffer(5:)
	 close(funit1)
 

	 funit1=findunitnum()
	 open(unit=funit1, file=avsfile, action='read')
	 read(unit=funit1,fmt=*) ncol 
	 isat=0 
	 ipor=0 
	 iden=0 
	 do i=1, ncol 
	  read(unit=funit1,fmt='(a)') buffer   
	  buffer=adjustl(buffer)
	  if( buffer(1:3) .eq. "Sat") isat=i  
	  if( buffer(1:3) .eq. "Por") ipor=i 
	  if( buffer(1:9) .eq. "Liquid De") iden=i 
	 end do
	 if( isat .eq. 0 .or. ipor .eq. 0 .or. iden .eq. 0) goto 996  
	 isat=isat+1 
	 ipor=ipor+1
	 iden=iden+1 

	 allocate(tmparr(ncol+1) ) 
	 nvals = ptr2mesh%ncells 
	 print *, nvals
	 do i=1,nvals 
	   read(unit=funit1, fmt=*) tmparr 
	   aobj%porosity(i)=tmparr(ipor)
	   aobj%saturation(i)=tmparr(isat) 
	   aobj%density(i)=tmparr(iden) 
	 end do 
	 deallocate(tmparr) 
	close(funit1) 
 endif

return 

999 print *, 'no control file found' 
    stop  
998 print  *, 'no dispersion tensor block' 
    stop 
997 print *, 'no avs file found' 
    stop  
996 print *, 'porosity,density or saturation not in avs file' 
    stop 
995 print *, 'error in zone file' 
    stop 

end function flowfield_  


subroutine reconstruct_velfield(ffobj)
	type(flowfield), target :: ffobj  
	type(mesh), pointer :: crntmesh 
	type(cell), pointer :: ptr2cell, ptr2nbr 
	real(kind=dp), dimension(:,:), pointer :: vfield  
	real(kind=dp), dimension(:,:,:), pointer :: disp
	real(kind=dp), dimension(:,:,:), pointer :: dcoef  
	real(kind=dp) :: absvel , lbar , alphaT, Dm  
	real(kind=dp), dimension(3) :: vel , avec 
	real(kind=dp), dimension(:,:), allocatable :: nrmls 
	real(kind=dp), dimension(:), allocatable :: fluxes  
	real(kind=dp), dimension(:), pointer :: dtdisp  
	integer :: ncells ,i,nn,m,j,n 
	real(kind=dp) :: small=1.0d-9 
	real(kind=dp) :: psd 


	crntmesh => ffobj%ptr2mesh 

	ncells = crntmesh%ncells

	vfield => ffobj%vfield
	disp => ffobj%disp 
	dcoef => ffobj%dcoef 
	dtdisp => ffobj%dtdisp 
	
	do i=1,ncells
		! first gather necessary info 

		ptr2cell => crntmesh%cells(i)

		nn=ptr2cell%nn
		nn=nn-1 ! self connection not included below

		allocate( nrmls(nn,3) )
		allocate( fluxes(nn) )
		m = 0
		do j=1,nn+1
			n = ptr2cell%nbrs(j)
			if( n .eq. i ) cycle ! ignore self connection
			m = m + 1
			ptr2nbr => crntmesh%cells(n)
			avec = ptr2cell%posi - ptr2nbr%posi
			 ! normalization not needed because
			 ! area array is actually area/distance
			 !d2 = sum(avec*avec)
			 !avec = avec/sqrt(d2)
			avec = avec * crntmesh%areas( ptr2cell%indx2area(j) )
			nrmls(m,:) = avec

			fluxes(m) = ffobj%fluxes( ptr2cell%indx2conn(j) )

		end do

		!  vel=vrecon( nrmls, fluxes, ptr2cell)  

		vel=vreconbound( nrmls, fluxes, ptr2cell)  


		psd=ffobj%porosity(i) * ffobj%saturation(i) * ffobj%density(i) 
		vel=vel/psd ! convert from flux to velocity  

		vfield(i,:)=vel
		
		deallocate( nrmls )
		deallocate( fluxes )
	end do
	
	
	

end subroutine reconstruct_velfield



!******************************************************
subroutine reconstruct_flowfield(ffobj) 

	type(flowfield), target :: ffobj  
	type(mesh), pointer :: crntmesh 
	type(cell), pointer :: ptr2cell, ptr2nbr 
	real(kind=dp), dimension(:,:), pointer :: vfield  
	real(kind=dp), dimension(:,:,:), pointer :: disp
	real(kind=dp), dimension(:,:,:), pointer :: dcoef  
	real(kind=dp) :: absvel , lbar , alphaT, Dm  
	real(kind=dp), dimension(3) :: vel , avec 
	real(kind=dp), dimension(:,:), allocatable :: nrmls 
	real(kind=dp), dimension(:), allocatable :: fluxes  
	real(kind=dp), dimension(:), pointer :: dtdisp  
	integer :: ncells ,i,nn,m,j,n 
	real(kind=dp) :: small=1.0d-9 
	real(kind=dp) :: psd 


	crntmesh => ffobj%ptr2mesh 	!get the mesh

	ncells = crntmesh%ncells 	!and the cells that the mesh defines

	vfield => ffobj%vfield		!velocity field
	disp => ffobj%disp 			!displacement matrix
	dcoef => ffobj%dcoef 		!dcoef matrix
	dtdisp => ffobj%dtdisp 		!time step limit by dispersive 

	if (.not.ffobj%amaread) then
		call reconstruct_velfield(ffobj)	!temporary call
	endif

	do i=1,ncells

		

		! dispersion coefficient
		disp(i,:,:)  = calculateBtensor(ffobj%tensor(i), vfield(i,:)) 
		dcoef(i,:,:) = calculateDtensor(ffobj%tensor(i), vfield(i,:)) 

		! local time step limits  
		absvel= sqrt( dot_product(vfield(i,:), vfield(i,:)) ) + small  
		call MaxTransverseDispersion( ffobj%tensor(i), i, alphaT, Dm) 

		lbar = crntmesh%vols(i)**(1.0d0/3.0d0) 
		dtdisp(i) = lbar*lbar/(4.0d0 * ( alphaT * absvel + Dm + small) )

		

	end do


	!  eventually we will have flux as a component of each fehm solution
	!  this can be deallocated once velocity field is built
	if (.not. ffobj%amaread)then
		deallocate( ffobj%fluxes )
	endif

end subroutine reconstruct_flowfield 

!*********************************************************************

function vrecon( amat, fluxes,ptr2cell ) result(vel)

  real(kind=dp), dimension(:,:) :: amat  
  real(kind=dp), dimension(:) :: fluxes
  type(cell), pointer :: ptr2cell 
  real(kind=dp), dimension(3) :: vel
  real(kind=dp), dimension(3,3) :: ata , atainv
  real(kind=dp), dimension(:,:), allocatable  ::  aaa
  integer :: m,n,nn

   nn=ptr2cell%nn-1

   ! now find velocity
   do m=1,3
     do n=1,3
       ata(m,n) = dot_product( amat(:,m), amat(:,n) )
     end do
   end do

   atainv=invert3x3(ata)

   allocate( aaa(3,nn) )
   do m=1,3
    do n=1,nn
      aaa( m, n) = dot_product( atainv(m,:),  amat(n,:) )
    end do
   end do

   do m=1,3
     vel(m) = dot_product( aaa(m,:), fluxes )
   end do
  
  deallocate(aaa) 

 return

end function vrecon

!*********************************************************************


function vreconbound( amat, fluxes,ptr2cell ) result(qhatb)

  real(kind=dp), dimension(:,:) :: amat
  real(kind=dp), dimension(:) :: fluxes
  type(cell), pointer :: ptr2cell
  real(kind=dp), dimension(3) :: qhat , qhatb 
  real(kind=dp), dimension(3,3) :: ata , atainv
  real(kind=dp), dimension(:,:), allocatable  ::  aaa
  real(kind=dp), dimension(:,:), allocatable  ::  atainvbt , batainvbt , btw , w 
  real(kind=dp), dimension(:,:), allocatable :: b 
  real(kind=dp), dimension(:), allocatable :: bqhat 
  integer :: m,n,nn,nb,nout,iout  
  logical :: inflow 

   nn=ptr2cell%nn-1

   ! gtranspose g  
   do m=1,3
     do n=1,3
       ata(m,n) = dot_product( amat(:,m), amat(:,n) )
     end do
   end do

   ! inverse(transpose(g) g) 
   atainv=invert3x3(ata)

   allocate( aaa(3,nn) )
   do m=1,3
    do n=1,nn
      aaa( m, n) = dot_product( atainv(m,:),  amat(n,:) )
    end do
   end do

   do m=1,3
     qhat(m) = dot_product( aaa(m,:), fluxes )
   end do

   ! finished with unconstrained reconstruction 

   nb=ptr2cell%NumNoflow 

   !done if not a boundary 
   if( nb .eq. 0 ) then 
       qhatb=qhat 
       return 
   end if 


   ! apply constraints only to outflow boundaries if inflow is not allowed 

   !V1.1 
   inflow=ptr2cell%inflow 

   nout=0 
   do m=1,nb 
      ! V1.1 
      if( .not. inflow) then 
         nout=nout+1 
      else if( dot_product( ptr2cell%bmat(m,:), qhat ) .gt. 0.0d0 ) then  
         nout=nout+1 
       end if 
   end do 

   if( nout .eq. 0) then 
       qhatb=qhat 
       return 
   end if 

   ! done if a corner 
   if( nout .eq. 3) then 
       qhatb=0.0d0 
       return 
   end if 

   
   allocate(b(nout,3) ) 
   iout=0 
   do m=1,nb 
      ! V1.1 
      if( .not. inflow) then 
         iout=iout+1 
         b(iout,:)=ptr2cell%bmat(m,:) 
      else if( dot_product( ptr2cell%bmat(m,:), qhat ) .gt. 0.0d0 ) then  
         iout=iout+1 
         b(iout,:)=ptr2cell%bmat(m,:) 
      end if 
   end do 


   nb=nout 

   allocate( bqhat(nb) ) 
   do m=1,nb 
      bqhat(m) = dot_product( b(m,:), qhat ) 
   end do 

   allocate( atainvbt(3,nb) ) 
   do m=1,3 
    do n=1,nb 
      atainvbt(m,n) = dot_product( atainv(m,:), b(n,:) ) 
    end do 
   end do 
  

   allocate( batainvbt(nb,nb) ) 
   do m=1,nb 
    do n=1,nb 
      batainvbt(m,n) = dot_product( b(m,:), atainvbt(:,n) ) 
    end do 
   end do 


   !invert batainvbt  
   allocate(w(nb,nb) ) 
   if( nb .eq. 1) w(1,1)=1.0d0/batainvbt(1,1) 
   if( nb .eq. 2) w=invert2x2( batainvbt ) 
 
   
   allocate(btw(3,nb) )     
   do m=1,3 
    do n=1,nb 
     btw(m,n) = dot_product( b(:,m), w(:,n) ) 
    end do 
   end do 

   deallocate( w) 
   allocate( w(3,nb) )      
   do m=1,3 
     do n=1,nb 
        w(m,n) = dot_product( atainv(m,:) , btw(:,n) ) 
     end do 
   end do 


   do m=1,3 
    qhatb(m) = qhat(m)-  dot_product( w(m,:) , bqhat ) 
   end do 


   deallocate(atainvbt) 
   deallocate( batainvbt) 
   deallocate(aaa) 
   deallocate(btw) 
   deallocate( w) 
   deallocate(bqhat)  
    

 return

end function vreconbound 


end module flowfield_class 

