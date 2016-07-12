program pclite 

 implicit none 

 character(len=30) :: buff 
 integer :: npart,part, idum  , icount  , ipart , icell, i
 integer, parameter :: dp=kind(1.0d0) 
 real(kind=dp) :: dum , tin, tout, treport  

 type phistory 
  real(kind=dp), dimension(:), pointer :: times 
  integer, dimension(:), pointer :: cells
 end type phistory 

 integer, dimension(:), pointer :: counts 
 integer :: ncells

 integer :: numargs 

 real(kind=dp), dimension(:), pointer :: score 
 real(kind=dp), dimension(:), pointer :: xpos, ypos, zpos 
 real(kind=dp), dimension(:), pointer :: vol  

 type(phistory), dimension(:), pointer :: tracks 

! get names file from command line
numargs = iargc()
if(numargs .eq. 0) then
 print *, 'please enter time at which concentration is required'  
 stop 
else
 call getarg(1, buff) 
 read(buff,*) treport  
end if

 open(unit=10, file='stor') 
 read(10,*) 
 read(10,*) 
 read(10,*) idum, ncells , idum, idum, idum 
 allocate( vol(ncells) ) 
 read(10,*) vol 
 close(10) 

 open(unit=10, file='fehmn') 
 read(10,*) 
 read(10,*) ncells 
 allocate(xpos(ncells)) 
 allocate(ypos(ncells)) 
 allocate(zpos(ncells)) 
 do icell=1,ncells 
  read(10,*) i, xpos(i),ypos(i),zpos(i) 
 end do 
 close(10) 
 
 open(unit=10, file='sptr2') 
 read(10,'(a)') 
 read(10,'(a)') 
 read(10,*) npart 
 read(10,'(a)') 

 allocate( tracks(npart) ) 
 allocate( counts(npart) ) 
 counts=0 

 do 
  read(unit=10, fmt=*, end=99) part, dum, idum  
  counts( part ) = counts(part)+1 
 end do 
99 continue 

 
 do ipart=1,npart 
   allocate( tracks(ipart)%times(0:counts(ipart)) )
   tracks(ipart)%times(0) = 0.0d0 
   allocate( tracks(ipart)%cells(counts(ipart)) ) 
 end do 


rewind(10) 
 read(10,'(a)') 
 read(10,'(a)') 
 read(10,*) npart 
 read(10,'(a)') 

 counts = 0 
 do 
  read(unit=10, fmt=*, end=999) ipart, dum, idum  
  if( ipart .gt. npart) exit ! allows analysis of subset by changing npart 
  counts(ipart) = counts(ipart) + 1 
  icount=counts(ipart) 
  tracks(ipart)%times(icount) = dum 
  tracks(ipart)%cells(icount) = idum 
 end do 
999 continue 
close(10) 


allocate(score(ncells)) 
score=0.0d0  
do ipart=1,npart 
 do icount=1,counts(ipart) 
   icell=tracks(ipart)%cells(icount) 
   tin=tracks(ipart)%times(icount-1) 
   tout=tracks(ipart)%times(icount) 
   tin=max(0.0, treport-tin) 
   tout=max(0.0, treport-tout) 
   score(icell) = score(icell) + tin-tout 
 end do 
end do 

do icell=1,ncells 
  print *, icell,xpos(icell),ypos(icell),zpos(icell),score(icell)/(float(npart) *vol(icell) ) 
end do 

end  program pclite 

