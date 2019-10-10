module mesh_class 

use precision  
use iotools 
use cell_class 
implicit none 

private 

public mesh,mesh_ 
public element 

type element 
 integer :: i,j,k,l 
 integer, dimension(4) :: nbrs 
end type element 

type mesh 
 integer :: ncells, nelems, nconns, nareas 
 type(cell), dimension(:), pointer :: cells 
 type(element), dimension(:), pointer :: elems 
 real(kind=dp),dimension(:), pointer :: areas 
 real(kind=dp),dimension(:), pointer :: vols  
end type mesh 

contains 

! mesh constructor *************************************
function mesh_( namefile )  result(themesh) 

character(len=70) :: namefile 
character(len=70) :: fehmnfile, storfile,eafile,bfile
character(len=70) :: ozfile
character(len=80) :: buffer 
type(mesh), target :: themesh  
integer :: count , icell 
integer :: funit 

 
 funit=findunitnum() 
 open(unit=funit, file=namefile,action='read') 

 do 
  read(unit=funit,fmt='(a)', end=999) buffer 
  buffer=adjustl(buffer) 
  buffer=trim(buffer) 
  if( buffer(1:6) .eq. "fehmn:") exit 
 end do 
 fehmnfile = buffer(7:) 
 call readfehmnfile(fehmnfile, themesh) 

 rewind(funit) 
 do 
  read(unit=funit,fmt='(a)', end=998) buffer 
  buffer=adjustl(buffer) 
  buffer=trim(buffer) 
  if( buffer(1:5) .eq. "stor:") exit 
 end do 
 storfile = buffer(6:) 
 call readstorfile(storfile, themesh) 

 rewind(funit) 
 do 
  read(unit=funit,fmt='(a)', end=997) buffer 
  buffer=adjustl(buffer) 
  buffer=trim(buffer) 
  if( buffer(1:7) .eq. "ealist:") exit 
 end do 
 eafile = buffer(8:) 
 call readeafile(eafile, themesh) 


 rewind(funit) 
 do 
  read(unit=funit,fmt='(a)', end=995) buffer 
  buffer=adjustl(buffer) 
  buffer=trim(buffer) 
  if( buffer(1:7) .eq. "cbound:") exit 
 end do 
 ozfile = buffer(8:) 
 print *, 'no-transport boundary file: ', ozfile 

 ! read the outside zone file 
 call readoutsidezonefile(ozfile,  themesh) 
 

 close(funit) 
 return 

999 print *, 'fehmn filename not found in name file'  
 stop 

998 print *, 'stor filename not found in name file'  
 stop 

997 print *, 'ealist filename not found in name file'  
 stop 

995 print *, 'no-flow boundaries are not defined'  
 do icell=1,themesh%ncells 
   themesh%cells(icell)%NumNoflow = 0 
 end do 
 close(funit) 
 return 

end function mesh_ 
! ******************************************************

! readeafile ****************************************
subroutine readeafile( eafile, themesh )  

type(mesh), target :: themesh  
character(len=70) :: eafile
type(element), pointer :: ptr1
integer :: ielem,nelems , i,j,k,l, idum, jdum  
integer :: funit 

funit=findunitnum() 
open(unit=funit, file=eafile, action='read') 

nelems=themesh%nelems 


do ielem=1,nelems 
 
 read(funit,*) idum, jdum, i,j,k,l 

 if( ielem .ne. idum .or. jdum .ne. 4) then 
  print *, 'problem reading element adjacency lists' 
  stop 
 end if 

 ptr1 => themesh%elems(ielem) 
 ptr1%nbrs(1:4) = (/i,j,k,l/)  

end do 
close(funit) 

end subroutine readeafile

! readfehmnfile ****************************************
subroutine readfehmnfile( fehmnfile, themesh )  

type(mesh), target :: themesh  
character(len=70) :: fehmnfile
type(cell), pointer :: p2c
type(element), pointer :: p2e
integer :: i,j,ncells,nelems  
real(kind=dp), dimension(3) :: posi 
integer :: funit 

character(len=30) :: key 


 ! get positions from fehmn file
 ! later to replace this with subroutines dynamically bound 
 !  for formatted and unformatted read 

 print *, 'opening file ', fehmnfile 
 funit=findunitnum() 
 open(unit=funit,file=fehmnfile,action='read')
 do
   read(unit=funit,fmt=*,end=666) key
   if(key .eq. 'coor') exit
 end do

 read(unit=funit,fmt=*) ncells 

 themesh%ncells=ncells 

 allocate( themesh%cells(ncells) )

 do i=1,ncells 
  read(unit=funit,fmt=*) j,posi
  if( i .ne. j ) goto 666  
  p2c => themesh%cells(i) 
  p2c%posi = posi
 end do

 rewind(funit) 

 ! now get elements 
 open(unit=funit,file=fehmnfile) 
 do
   read(unit=funit,fmt=*,end=666) key
   if(key .eq. 'elem') exit
 end do

 read(unit=funit,fmt=*) j,nelems 
 if(j .ne. 4) goto 666 
 
 themesh%nelems = nelems 
 allocate( themesh%elems(nelems) ) 
 do i=1,nelems 
  p2e => themesh%elems(i) 
  read(unit=funit,fmt=*) j,p2e%i, p2e%j, p2e%k, p2e%l 
  if(i .ne. j) goto 666 
 end do 

 close(funit)

 return 

666 continue 
  print *, 'error in mesh constructor' 
  stop 

end subroutine readfehmnfile   

!*******************************************************************

!readstorfile ******************************************************
subroutine readstorfile( storfile, themesh )  

! reads store file produced by LaGrit 
! upon entry some components of themesh have already been allocated
!   including cells, elements
! areas have not been allocated 

type(mesh), target :: themesh  
character(len=70) :: storfile
type(cell), pointer :: p2c
integer, dimension(:), allocatable :: ncarr ,parray,idum ,  indxarea
integer :: funit 

integer :: i,j,m,n ,nn,nareas 

integer :: ncells,nsize,nconns 

 ! get connections from stor files
 funit=findunitnum() 
 print *, 'opening file ' , storfile 
 open(unit=funit,file=storfile, action='read')

 read(unit=funit,fmt=*)
 read(unit=funit,fmt=*)
 read(unit=funit,fmt=*) nareas, ncells, nsize

 if( ncells .ne. themesh%ncells ) goto 666  
 allocate(themesh%vols(ncells))
 read(unit=funit,fmt=*) themesh%vols 

 nconns=nsize-ncells-1
 themesh%nconns = nconns 
 themesh%nareas = nareas 
 allocate(ncarr(nconns) )
 allocate(parray(ncells+1) )
 read(unit=funit,fmt=*) parray
 read(unit=funit,fmt=*) ncarr

 allocate( indxarea(nsize) )   ! is padded by nnodes+1. will not use last values 
 read(unit=funit, fmt=*) indxarea

 allocate( idum(ncells) ) ! location of diagonal in file, but not needed
 read(unit=funit, fmt=*) idum
 deallocate( idum)

 allocate(themesh%areas(nareas) )
 read(unit=funit, fmt=*) themesh%areas 

 close(funit)
 ! finished with read

 ! number neighbors for each cell  
 themesh%ncells = ncells 
 do i=1,ncells 
  p2c => themesh%cells(i) 
  p2c%nn=parray(i+1)-parray(i)
 end do

 ! load the neighbor list
 m=0
 do i=1,ncells 
  p2c => themesh%cells(i) 
  nn = p2c%nn 
  n=0
  allocate( p2c%nbrs(nn) )  
  allocate( p2c%indx2area(nn) )  
  allocate( p2c%indx2conn(nn) )  
  do j=1,nn
    m=m+1
    n=n+1
    p2c%nbrs(n) = ncarr(m) 
    p2c%indx2area(n) = indxarea(m) 
    p2c%indx2conn(n) = m !check on this 
  end do
 end do

 deallocate(indxarea) 
 deallocate(parray)
 deallocate(ncarr)

 return 

666 continue 

 print *, ' error in readstorfile'  
 stop 

end subroutine readstorfile 

!********************************************************
 function isneighbor(a,b) result( isnbr )
 ! dead code, not used  

 type(element), pointer :: a, b
 integer, dimension(4) :: avec, bvec
 integer :: i,j, shared 
 logical :: isnbr

    avec=(/ a%i,a%j,a%k,a%l /)
    bvec=(/ b%i,b%j,b%k,b%l /)

    isnbr=.false. 
    shared=0
    do i=1,4
     do j=1,4
      if( avec(i) .eq. bvec(j) ) then 
        shared=shared+1
        exit 
      end if 
     end do
     if( i.eq.2 .and. shared .lt. 1) exit 
     if( i.eq.3 .and. shared .lt. 2) exit 
    end do
    if( shared .eq. 3 ) isnbr=.true.

 end function isneighbor 

!*********************************************************************
subroutine readoutsidezonefile( zonefile, themesh) 

type(mesh), target :: themesh 
character(len=70) :: zonefile,  buffer 
character(len=7), dimension(30) :: facestring 
integer ::  nnum , iface 
integer :: icell, ielem, i, j, k, l  
logical :: errflag 
integer, allocatable, dimension(:) :: blist  
integer, allocatable, dimension(:) :: counts  
type(cell), pointer :: thecell 
real(kind=dp), dimension(3) :: avec 
type(cell), pointer :: ptr2cell 
type(cell), pointer :: cell1, cell2, cell3, cell4  
type(element),pointer :: ptr2elem
integer :: funit 
logical :: flag1, flag2, flag3, flag4 
logical :: inflow   

funit=findunitnum() 
open(unit=funit, file=zonefile, action='read') 

!V1.1 determine how to treat the boundaries 
!  currently, this is global but will be broadcast to the cells  
!  could make specific to each face in the future 

inflow=.true. 
read(unit=funit,fmt='(a)') buffer
if( index(buffer, 'NOFLOW') .ne. 0) inflow=.false. 


facestring(1)='top' 
facestring(2)='bottom' 
facestring(3)='left_w' 
facestring(4)='right_e' 
facestring(5)='back_n' 
facestring(6)='front_s' 

! first determine number of boundary faces for each cell 

allocate(counts(themesh%ncells) ) 
counts=0

do iface=1,6

  nnum=0
  do
   read(unit=funit,fmt='(a)',end=997) buffer
   if( index(buffer, adjustl(facestring(iface) )) .ne. 0) exit
  end do
  read(unit=funit, fmt='(a)') buffer
  read(unit=funit, fmt=*) nnum

  allocate( blist(nnum) ) 

  read(unit=funit, fmt=*) blist

  do i=1,nnum 
   counts( blist(i) ) = counts( blist(i) ) + 1 
  end do 

  deallocate(blist) 
 
997 rewind(funit)

end do  ! end loop over face types
close(funit) 

! at this point we know the number of boundary faces for each cell 
! allocate and load the bmats 

do icell=1,themesh%ncells 
  thecell => themesh%cells(icell) 
  thecell%NumNoflow= counts(icell) 
  thecell%inflow=inflow 
  allocate( thecell%bmat( counts(icell), 3 ) ) 
end do 


funit=findunitnum() 
open(unit=funit, file=zonefile, action='read')

counts=0 ! will reuse as running count  
do iface=1,6
  
  nnum=0
  do
   read(unit=funit,fmt='(a)',end=996) buffer
   if( index(buffer, adjustl(facestring(iface) )) .ne. 0) exit
  end do
  read(unit=funit, fmt='(a)') buffer
  read(unit=funit, fmt=*) nnum

  allocate( blist(nnum) )

  read(unit=funit, fmt=*) blist

  do i=1,nnum 

   icell = blist(i) 
   counts(icell) = counts(icell) + 1    
   select case(iface) 
     case(1) 
      avec=(/0,0,1/) 
     case(2) 
      avec=(/0,0,-1/) 
     case(3)   
      avec=(/-1,0,0/) 
     case(4)  
      avec=(/1,0,0/)   
     case(5)  
      avec=(/0,1,0/) 
     case(6)  
      avec=(/0,-1,0/) 
   end select 
   themesh%cells(icell)%bmat(counts(icell), :) = avec 

  end do 

  deallocate(blist) 
 
996 rewind(funit)

end do  ! end loop over face types
close(funit) 

deallocate(counts) 


! set the no-transport boundaries 
! boundaries that are no-flow are given the designator 0 in the 
!   element neighbor list 
! boundaries that are open are given the designator -99  

  do ielem=1, themesh%nelems 

        ptr2elem => themesh%elems(ielem) 

        i = ptr2elem%i
        j = ptr2elem%j
        k = ptr2elem%k
        l = ptr2elem%l

        flag1=.false. 
        flag2=.false. 
        flag3=.false. 
        flag4=.false. 
        if( themesh%cells(i)%NumNoflow .gt. 0) flag1=.true.
        if( themesh%cells(j)%NumNoflow .gt. 0) flag2=.true.
        if( themesh%cells(l)%NumNoflow .gt. 0) flag3=.true. !note reversal of l and k 
        if( themesh%cells(k)%NumNoflow .gt. 0) flag4=.true.

        if( ptr2elem%nbrs(1) .le. 0) ptr2elem%nbrs(1) = -99
        if( ptr2elem%nbrs(2) .le. 0) ptr2elem%nbrs(2) = -99
        if( ptr2elem%nbrs(3) .le. 0) ptr2elem%nbrs(3) = -99
        if( ptr2elem%nbrs(4) .le. 0) ptr2elem%nbrs(4) = -99
!        if( flag2 .and. flag3 .and. flag4 .and. ptr2elem%nbrs(1) .lt. 0) ptr2elem%nbrs(1) = 0
!        if( flag1 .and. flag4 .and. flag3 .and. ptr2elem%nbrs(2) .lt. 0) ptr2elem%nbrs(2) = 0
!        if( flag1 .and. flag2 .and. flag4 .and. ptr2elem%nbrs(3) .lt. 0) ptr2elem%nbrs(3) = 0
!        if( flag1 .and. flag3 .and. flag2 .and. ptr2elem%nbrs(4) .lt. 0) ptr2elem%nbrs(4) = 0

        cell1 => themesh%cells(i) 
        cell2 => themesh%cells(j) 
        cell3 => themesh%cells(l) ! note reversal 
        cell4 => themesh%cells(k) 

        if( flag2 .and. flag3 .and. flag4 .and. ptr2elem%nbrs(1) .lt. 0) then 
            ! if here then nodes on facet are all on a NoFlow boundary, but do they share the same one? 
            if( shared_facetQ( cell2%bmat, cell3%bmat, cell4%bmat) )  ptr2elem%nbrs(1) = 0
        end if 
        if( flag1 .and. flag4 .and. flag3 .and. ptr2elem%nbrs(2) .lt. 0) then 
            if( shared_facetQ( cell1%bmat, cell4%bmat, cell3%bmat) )  ptr2elem%nbrs(2) = 0
        end if 
        if( flag1 .and. flag2 .and. flag4 .and. ptr2elem%nbrs(3) .lt. 0) then 
            if( shared_facetQ( cell1%bmat, cell2%bmat, cell4%bmat) )  ptr2elem%nbrs(3) = 0
        end if 
        if( flag1 .and. flag3 .and. flag2 .and. ptr2elem%nbrs(4) .lt. 0) then 
            if( shared_facetQ( cell1%bmat, cell3%bmat, cell2%bmat) )  ptr2elem%nbrs(4) = 0
        end if 
   
  end do 

return 
end subroutine readoutsidezonefile 

function shared_facetQ( amat, bmat, cmat) result(flag) 

 real(kind=dp), dimension(:,:) :: amat  , bmat, cmat 
 logical :: flag 
 integer :: n1, n2, n3, i, j, k, count   

 flag=.false. 

 n1=size(amat, 1) 
 n2=size(bmat, 1) 
 n3=size(cmat, 1) 
 
 count=0 
 do i=1,n1 
  do j=1,n2 
   if(abs( dot_product( amat(i,:), bmat(j,:) ) -1.0d0) .lt. 1.0d-4) then  
     do k=1,n3 
      if(abs( dot_product( amat(i,:), cmat(k,:) ) -1.0d0) .lt. 1.0d-4) count=count+1  
     end do 
   end if 
  end do 
 end do 
!  if( count .gt. 1) then 
!    print *, 'error in shared_facetQ' 
!    print *, count 
!    print *, n1
!    print *, amat 
!    print *, n2 
!    print *, bmat 
!    print *, n3 
!    print *, cmat 
!    stop 
!  end if 
!  if( count .eq. 1) flag=.true. 
  if( count .ge. 1) flag=.true. 

end function  

end module mesh_class 

