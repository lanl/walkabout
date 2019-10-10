module tracktools 

      use precision 
      use cell_class 
      use mesh_class 
      use linked_list  
      implicit none

private 

public in_elem 
public interpolate_tet1d   ! redundant capability tet_d is more general 
public interpolate_tet3d 
public invert3x3 
public invert2x2 
public findelement 
public findcell 
public intersect_face 
public initsearch 
public spiralsearch 

integer, dimension(:), allocatable :: visited  ! fix time permits.  Will not work if using multiple meshes 

contains 

!**********************************************************************
function findcell( crntmesh, ptr2elem, pt) result(icell) 

! finds vertex in the tet-element closest to pt 

integer :: icell 
integer :: i, j, k, l , inbr 
real(kind=dp), dimension(3) :: vert1, vert2, vert3, vert4 ,  pt , dpt
real(kind=dp) :: dmin, d2 
type(element), pointer :: ptr2elem , ptr2nbr 
type(mesh), pointer :: crntmesh 

   dmin=1.0d20     

   do inbr=0,4 

    if( inbr .eq. 0) then  
      ptr2nbr => ptr2elem 
    else 
      if( ptr2elem%nbrs(inbr)  .le. 0) cycle 
      ptr2nbr => crntmesh%elems( ptr2elem%nbrs(inbr) ) 
    end if 

    i=ptr2nbr%i
    j=ptr2nbr%j
    k=ptr2nbr%k
    l=ptr2nbr%l
    vert1 = crntmesh%cells( i )%posi
    vert2 = crntmesh%cells( j )%posi
    vert3 = crntmesh%cells( k )%posi
    vert4 = crntmesh%cells( l )%posi

   dpt = pt-vert1 
   d2 = dot_product( dpt, dpt) 
   if( d2 .lt. dmin) then 
     dmin=d2 
     icell=i 
   end if  

   dpt = pt-vert2 
   d2 = dot_product( dpt, dpt) 
   if( d2 .lt. dmin) then 
     dmin=d2 
     icell=j 
   end if 
  
   dpt = pt-vert3 
   d2 = dot_product( dpt, dpt) 
   if( d2 .lt. dmin) then 
     dmin=d2 
     icell=k 
   end if 

   dpt = pt-vert4 
   d2 = dot_product( dpt, dpt) 
   if( d2 .lt. dmin) then 
     dmin=d2 
     icell=l 
   end if 

   end do 


end function findcell 

! ######################################################################

   function findelement( crntmesh,  pnt) result(ielem) 

   type(mesh), pointer :: crntmesh 
   real(kind=dp), dimension(3) :: pnt 
   integer :: ielem,nelems 
   logical :: iflag 

   nelems=crntmesh%nelems 

   do ielem=1,nelems 
   
    iflag=in_elem( crntmesh, ielem, pnt) 
    if(iflag) return

   end do 

   ! if here the was not found 
   ielem=-1 
   return 

   end function findelement 

! ######################################################################
   
   function in_elem(crntmesh, m, pnt) result(isin) 
   type(mesh), pointer :: crntmesh 
   logical :: isin 
   integer :: m, iflag  
   real(kind=dp), dimension(3) :: pnt 

      iflag=inside_element(crntmesh,m,pnt) 
      isin=.false. 
      if(iflag .ge. 0) isin=.true. 

   end function in_elem 

! #####################################################################
   function inside_element( crntmesh, m, pnt) result(iflag) 



   type(mesh), pointer :: crntmesh 
   integer :: m 

   type(cell), pointer :: crntcell 

   real(kind=dp), dimension(3) :: pnt 
   real(kind=dp), dimension(:),  pointer :: pnt1, pnt2, pnt3, pnt4

   integer :: iflag,  i, j, k, l 

 
      i= crntmesh%elems(m)%i  
      j= crntmesh%elems(m)%j  
      k= crntmesh%elems(m)%k  
      l= crntmesh%elems(m)%l  

      crntcell => crntmesh%cells(i) 
      pnt1 => crntcell%posi 

      crntcell => crntmesh%cells(j) 
      pnt2 => crntcell%posi 

      crntcell => crntmesh%cells(l) !note reversal of l and k 
      pnt3 => crntcell%posi 

      crntcell => crntmesh%cells(k) 
      pnt4 => crntcell%posi 


      iflag=inside_tet(pnt1, pnt2, pnt3, pnt4, pnt) 

   end function inside_element 

   function inside_tet(pnt1, pnt2, pnt3, pnt4, pnt) result(iflag) 
   real(kind=dp), dimension(3) :: pnt , pnt1, pnt2, pnt3, pnt4
   integer ::  iflag234,iflag143,iflag124,iflag132, iflag 
   real(kind=dp) :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
   real(kind=dp) :: xa,ya,za

   ! Scott Painter, modified from LaGrit routines provided by Carl Gable 
  
      x1 = pnt1(1)  
      y1 = pnt1(2) 
      z1 = pnt1(3)  

      x2 = pnt2(1)  
      y2 = pnt2(2) 
      z2 = pnt2(3)  

      x3 = pnt3(1)  
      y3 = pnt3(2) 
      z3 = pnt3(3)  

      x4 = pnt4(1)  
      y4 = pnt4(2) 
      z4 = pnt4(3)  

      xa = pnt(1)  
      ya = pnt(2) 
      za = pnt(3)  

      iflag234=0
      iflag143=0
      iflag124=0
      iflag132=0
      call inside_tri(x2,y2,z2,x3,y3,z3,x4,y4,z4,xa,ya,za,iflag234)
      call inside_tri(x1,y1,z1,x4,y4,z4,x3,y3,z3,xa,ya,za,iflag143)
      call inside_tri(x1,y1,z1,x2,y2,z2,x4,y4,z4,xa,ya,za,iflag124)
      call inside_tri(x1,y1,z1,x3,y3,z3,x2,y2,z2,xa,ya,za,iflag132)
      if(iflag234.eq.0.and.iflag143.eq.0 .and.  iflag124.eq.0.and.iflag132.eq.0) then
         iflag=0  ! inside element 
      elseif(iflag234.ge.0.and.iflag143.ge.0 .and. iflag124.ge.0.and.iflag132.ge.0) then  ! on a face 
         if(iflag234.gt.0) then
            iflag=1  ! on face 1 
         elseif(iflag143.gt.0) then
            iflag=2  ! on face 2 
         elseif(iflag124.gt.0) then
            iflag=3  ! on face 3 
         elseif(iflag132.gt.0) then
            iflag=4  ! on face 4 
         else
            iflag=-1
         endif
      else
         iflag=-1   ! point is outside element 
      endif

      return
      end function inside_tet 

      subroutine inside_tri(xl1,yl1,zl1,xl2,yl2,zl2,xl3,yl3,zl3,  & 
                           xa,ya,za,iflag)
! routine comes from LaGrit system 
! ######################################################################
!
      implicit none
! 
! ######################################################################
!
      real(kind=dp) xl1,yl1,zl1, xl2,yl2,zl2, xl3,yl3,zl3
      real(kind=dp) xa,ya,za
      real(kind=dp) ax4,ay4,az4
      real(kind=dp)  xfacbox
      real(kind=dp) xminbox, yminbox, zminbox
      real(kind=dp) xmaxbox, ymaxbox, zmaxbox
      integer iflag

      real(kind=dp) voltet

      real(kind=dp) xepsilon
      data xepsilon / 1.0d-07 /
!
!  ######################################################################
!
      xminbox=min(xl1,xl2,xl3)
      yminbox=min(yl1,yl2,yl3)
      zminbox=min(zl1,zl2,zl3)
      xmaxbox=max(xl1,xl2,xl3)
      ymaxbox=max(yl1,yl2,yl3)
      zmaxbox=max(zl1,zl2,zl3)
      xfacbox=xepsilon*sqrt((xmaxbox-xminbox)**2 +  & 
                       (ymaxbox-yminbox)**2 +      & 
                       (zmaxbox-zminbox)**2 )
      AX4=  (YL3-YL1)*(ZL2-ZL1)-(ZL3-ZL1)*(YL2-YL1)
      AY4=-((XL3-XL1)*(ZL2-ZL1)-(ZL3-ZL1)*(XL2-XL1))
      AZ4=  (XL3-XL1)*(YL2-YL1)-(YL3-YL1)*(XL2-XL1)
      VOLTET=-((xa-XL1)*AX4+(ya-YL1)*AY4+(za-ZL1)*AZ4)
!
!     ******************************************************************
!
!     IFLAG == +1 ==> POINT MUST BE ON SURFACE.
!           ==  0 ==> POINT MUST BE BELOW SURFACE.
!           == -1 ==> POINT MUST BE ABOVE SURFACE.
!
      if(abs(voltet).le.xfacbox) then
         iflag=+1
      elseif(voltet.gt. xfacbox) then
         iflag=-1
      elseif(voltet.lt.-xfacbox) then
         iflag= 0
      endif

      return
      end subroutine inside_tri

!*************************************************************************************
      function interpolate_tet1d(vert1, vert2, vert3, vert4, vertvals, pt, idir) result( aresult )  
      ! interpolates a scalar on a tet using barycentric coordinates  
      
      real(kind=dp), dimension(3) :: vert1, vert2, vert3, vert4 ,  pt , apt 
      real(kind=dp) :: aresult  
      real(kind=dp), dimension(3,3) :: amat , amatinv 
      real(kind=dp), dimension(4) :: lam
      real(kind=dp), dimension(4) :: vertvals  

      integer :: i , idir 

      amat(1,:)=vert1-vert4 
      amat(2,:)=vert2-vert4 
      amat(3,:)=vert3-vert4 

      amatinv=invert3x3(amat) 

      do i=1,3
       select case (idir)
        case(0)
          apt=pt - vert4
        case(1)
          apt=(/1,0,0/)
        case(2)
          apt=(/0,1,0/)
        case(3)
          apt=(/0,0,1/)
       end select
       lam(i) = dot_product( amatinv(:,i), apt)
      end do

      if( idir .eq. 0) then
      lam(4)=1.0-lam(1)-lam(2)-lam(3)
      else
      lam(4) = -lam(1)-lam(2)-lam(3)
      end if

      aresult =  dot_product( vertvals(:), lam) 
      
      return   

      end function interpolate_tet1d 

!*************************************************************************************
      function interpolate_tet3d(vert1, vert2, vert3, vert4, vertvels, pt,idir ) result( vel )

      !interpolates a 3-component vector on a tet 
      ! use for velocity 

      real(kind=dp), dimension(3) :: vert1, vert2, vert3, vert4 ,  pt
      real(kind=dp), dimension(3) :: vel, apt 
      real(kind=dp), dimension(3,3) :: amat , amatinv
      real(kind=dp), dimension(4) :: lam
      real(kind=dp), dimension(4,3) :: vertvels
    

      integer :: i,idir 

      amat(1,:)=vert1-vert4
      amat(2,:)=vert2-vert4
      amat(3,:)=vert3-vert4

      amatinv=invert3x3(amat)

      do i=1,3
       select case (idir)    
        case(0)  
          apt=pt - vert4 
        case(1) 
          apt=(/1,0,0/) 
        case(2) 
          apt=(/0,1,0/) 
        case(3)  
          apt=(/0,0,1/) 
       end select 
       lam(i) = dot_product( amatinv(:,i), apt)
      end do

      if( idir .eq. 0) then 
      lam(4)=1.0-lam(1)-lam(2)-lam(3)
      else 
      lam(4) = -lam(1)-lam(2)-lam(3) 
      end if 

!      if( ( any(lam .gt. 1.0d0) .or. any(lam .lt. 0.0) ) .and. idir .eq. 0) then
!         print *, 'warning: weight is less than zero ', lam
!         print *, vert1
!         print *, vert2
!         print *, vert3
!         print *, vert4
!         print *, pt
!      end if

      vel(1) =  dot_product( vertvels(:,1), lam)
      vel(2) =  dot_product( vertvels(:,2), lam)
      vel(3) =  dot_product( vertvels(:,3), lam)


      end function interpolate_tet3d 

!******************************************************************


   function invert2x2(a) result(ainv)

      real(kind=dp), dimension(2,2) :: a,ainv

      real(kind=dp) :: det

      det = a(1,1)*a(2,2)- a(1,2)*a(2,1) 

      ainv(1,1) =  a(2,2)  
      ainv(2,1) =  -a(2,1) 

      ainv(1,2) = -a(1,2)  
      ainv(2,2) = a(1,1)  

      ainv = ainv/det

   end function invert2x2

!***************************************************************


!******************************************************************


   function invert3x3(a) result(ainv)
   ! carl gable 

      real(kind=dp), dimension(3,3) :: a,ainv

      real(kind=dp) :: det

      det = a(1,1)*( a(3,3)*a(2,2) - a(3,2)*a(2,3) ) &
         &  - a(2,1)*( a(3,3)*a(1,2) - a(3,2)*a(1,3) ) &
         &  + a(3,1)*( a(2,3)*a(1,2) - a(2,2)*a(1,3) )

      ainv(1,1) =    a(3,3)*a(2,2) - a(3,2)*a(2,3)
      ainv(2,1) = -( a(3,3)*a(2,1) - a(3,1)*a(2,3) )
      ainv(3,1) =    a(3,2)*a(2,1) - a(3,1)*a(2,2)

      ainv(1,2) = -( a(3,3)*a(1,2) - a(3,2)*a(1,3) )
      ainv(2,2) =    a(3,3)*a(1,1) - a(3,1)*a(1,3)
      ainv(3,2) = -( a(3,2)*a(1,1) - a(3,1)*a(1,2) )

      ainv(1,3) =    a(2,3)*a(1,2) - a(2,2)*a(1,3)
      ainv(2,3) = -( a(2,3)*a(1,1) - a(2,1)*a(1,3) )
      ainv(3,3) =    a(2,2)*a(1,1) - a(2,1)*a(1,2)

      ainv = ainv/det

   end function invert3x3

!***************************************************************

   subroutine intersect_face ( crntmesh, m, pt1, pt2, bflag, zgb) 

   ! checks for intersection on face and returns 
   !  new element (m) on other side of face intersection and new point (pt1) 
   ! algorithm works recursively until m contains pt2 
   ! in the event that pt2 is outside the domain, pt1 is reset to the boundary 
   ! and m is set to the last element visited. flag is set if on boundary 

   type(mesh), pointer :: crntmesh 
   type(cell), pointer :: node1,node2, node3, node4

   real(kind=dp) :: d1, d2, d3, d4  ,dmin 
   real(kind=dp), dimension(3) :: pt1, pt2, pt  
   real(kind=dp), parameter :: small=1.0d-9  

   logical :: bflag  , inelem , zgb 

   integer :: iflag,  m , i, j, k, l , iface ,imin, nextelem 

      bflag=.false. ! boundary flag  
      zgb = .false. ! zero gradient boundary 

1001  continue 
         
      !v1.2 moved here from above 
      if( in_elem(crntmesh, m, pt2) ) return 

      i= crntmesh%elems(m)%i  
      j= crntmesh%elems(m)%j  
      k= crntmesh%elems(m)%k  
      l= crntmesh%elems(m)%l  

      node1 => crntmesh%cells(i) 
      node2 => crntmesh%cells(j) 
      node3 => crntmesh%cells(l) !note reversal of l and k. do not change! 
      node4 => crntmesh%cells(k) 

      iface = 0 

      call lineseg_tri(node2%posi, node3%posi, node4%posi, pt1, pt2,pt, iflag)
      if( iflag .eq. 0) then 
         iface=1 
         goto 998 
      end if  

      call lineseg_tri(node1%posi, node3%posi, node4%posi, pt1, pt2,pt, iflag)
      if( iflag .eq. 0) then 
         iface=2 
         goto 998 
      end if  

      call lineseg_tri(node1%posi, node2%posi, node4%posi, pt1, pt2,pt, iflag)
      if( iflag .eq. 0) then 
         iface=3 
         goto 998 
      end if  

      call lineseg_tri(node1%posi, node2%posi, node3%posi, pt1, pt2,pt, iflag)
      if( iflag .eq. 0) then 
         iface=4 
         goto 998 
      end if  


998 continue 

      nextelem=m  
      if( iface .ne. 0 ) then  
       nextelem=crntmesh%elems(m)%nbrs(iface) 
      end if 


      if( nextelem .ne. m .and. nextelem .gt. 0) then 
          pt1 = pt 
          m=nextelem       
          goto 1001 
      end if 

      if( nextelem .le. 0) then 
        bflag=.true. 
        pt1 = pt
	if( nextelem .eq. 0) zgb=.true. 
        return 
      end if 

! check to make sure we are where we think we should be 

      m=nextelem 
      inelem=in_elem(crntmesh, nextelem, pt2)
      if( .not. inelem) then 
        nextelem=spiralsearch(crntmesh, m, pt2) 
        if( nextelem .le. 0) then 
         bflag=.true. 
         zgb=.true. 
         print *, 'searching for cell ...' 
        else 
         m=nextelem 
        end if 
      end if 
	  

      return 

   end subroutine intersect_face



!***************************************************************

      function line_intersect_plane(x1,x2,pt1,pt2,pt3) result(tt) 
      ! find point of intersection of the line given by x1;x2
      ! and the plane given by nodes i1,i2,i3

      ! from LaGrit or maybe SPTR?

      implicit none

      integer i

      real(kind=dp), dimension(3) :: x1,x2, pt1, pt2, pt3 
      real(kind=dp) ::  cx,cy,cz,cc,cl,cm,cn
      real(kind=dp), dimension(3,3)  :: a(3,3), coor(3,3) 
      real(kind=dp) ::  epsilon,fac1,fac2,tt

      epsilon = 1.e-18

! equation of the plane, Solid Analytical Geom by John Olmsted, 
! call #: QA553.O45.c7.  x*cx+y*cy+z*cz+cc=0.

      coor(1,:)=pt1 
      coor(2,:)=pt2 
      coor(3,:)=pt3 


      a(1,3)=1.
      a(2,3)=1.
      a(3,3)=1.

      do i=1,3
         a(i,1)=coor(i,2)
         a(i,2)=coor(i,3)
      enddo
      cx=det3(a)

      do i=1,3
         a(i,1)=coor(i,1)
         a(i,2)=coor(i,3)
      enddo
      cy=det3(a)
      cy=-cy

      do i=1,3
         a(i,1)=coor(i,1)
         a(i,2)=coor(i,2)
      enddo
      cz=det3(a)

      do i=1,3
         a(i,1)=coor(i,1)
         a(i,2)=coor(i,2)
         a(i,3)=coor(i,3)
      enddo
      cc=det3(a)
      cc=-cc

! equations of the line thru x1-x2

      cl=x2(1)-x1(1)
      cm=x2(2)-x1(2)
      cn=x2(3)-x1(3)

! parametric solution

      fac1=cx*cl+cy*cm+cz*cn

      if(abs(fac1).lt.epsilon) then
! this means the line is parallel to the plane and they do not
! intersect
         tt=-1.e+20

      else
         fac2=-(cx*x1(1)+cy*x1(2)+cz*x1(3)+cc)
         tt=fac2/fac1
      endif

!      xc(1)=x1(1)+cl*tt
!      xc(2)=x1(2)+cm*tt
!      xc(3)=x1(3)+cn*tt
      
      return

end function line_intersect_plane

! ************************************************************

function det3(a) result(deta) 

      implicit none

      real(kind=dp) :: a(3,3),deta

      deta=0.
      deta=deta+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      deta=deta-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      deta=deta+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

return

end function det3 


!*********************************************************
function spiralsearch(crntmesh,oldelem,pt) result(melem)  
! breadth-first search, starting at oldelem 
! returns index of element containing pt 

integer :: oldelem,melem ,ii ,m , isearched 
logical :: inelem 
real(kind=dp), dimension(3) :: pt 
type(mesh), pointer :: crntmesh 

    if(.not.allocated(visited)) then 
     allocate(visited(crntmesh%nelems))
     visited=0
    end if 

    melem=oldelem 

    ! start spiral search
    inelem = .false.
    call add2queuend(melem)
    visited(melem) = 1
    do
       melem=fromqueue( )
       if( melem .eq. -1) exit

       inelem=in_elem(crntmesh, melem, pt)
       if( inelem ) exit  ! success
       do ii=1,4
         m=crntmesh%elems(melem)%nbrs(ii)
         if( m .le. 0) cycle
         if( visited(m) .eq. 0) call add2queuend(m)
         visited(m) = 1
       end do

    end do

    ! now clear the queue
    if( melem .ne. -1) visited(melem) = 0
    isearched=0 
    do
       ii=destroyqueue()
       if( ii .eq. -1) exit
       isearched=isearched+1
       visited(ii) = 0
    end do
!    print *, 'searched ', isearched 

end function spiralsearch

!******************************************************
subroutine initsearch(crntmesh) 
type(mesh), pointer :: crntmesh 

if(.not.allocated(visited)) allocate(visited(crntmesh%nelems)) 
visited=0 

end subroutine initsearch 

  subroutine lineseg_tri(pt1, pt2, pt3,  pta, ptb, pt, iflag) 
      implicit none
      real(kind=dp) local_eps,releps,xa,xb,ya,yb,za,zb,dist1,dist2,det, & 
       xsum,area31,ax31,ay31,az31,area23,ax23,ay23,az23,area12, & 
       area123,ax12,ay12,az12,ax123,ay123,az123,dsum,dsb,dsa,  & 
       xunit,yunit,zunit,dsab,b3,a31,a32,a33,b2,a21,a22,a23,   & 
       b1,a13,a12,a11,cx,cy,cz,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(kind=dp) a1norm, a2norm, a3norm, d1, d2 

      real(kind=dp), dimension(3) :: pt1, pt2 , pt3 ,  pta, ptb , pt 

      integer iflag, i, ierror
      real*8 xnodes(3),ynodes(3),znodes(3)
      real*8 lennorm
      real*8 triarea
      character*132 warnmess

      x1=pt1(1) 
      x2=pt2(1) 
      x3=pt3(1) 

      y1=pt1(2) 
      y2=pt2(2) 
      y3=pt3(2) 

      z1=pt1(3) 
      z2=pt2(3) 
      z3=pt3(3) 

      xa=pta(1) 
      xb=ptb(1) 

      ya=pta(2) 
      yb=ptb(2) 

      za=pta(3) 
      zb=ptb(3) 
!
!C ######################################################################
!
      local_eps = 1e-8
!
!
!     Calculate the plane defined by the triangle
      cx=  (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
      cy=-((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1))
      cz=  (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
      triarea=0.5*sqrt(cx**2+cy**2+cz**2)
!
!     Releps = local_eps*.5*(length of line+area of triangle)
      releps = local_eps*(sqrt((xa-xb)**2+(ya-yb)**2+(za-zb)**2)+  & 
              triarea)
 
!
!     Create a matrix with the first row the equation of the plane
!     of the triangle, the second and third rows the equations of
!     two planes that make up the line.
!
      a11=cx
      a12=cy
      a13=cz
!     Normalize the first row of the matrix to length 1
      a1norm=sqrt(a11**2+a12**2+a13**2)
      a11=a11/a1norm
      a12=a12/a1norm
      a13=a13/a1norm

!     Do the point on the plane
      b1=(cx*x1+cy*y1+cz*z1)/a1norm

      a21=0.0
      a22=0.0
      a23=0.0
      b2=0.0
      a31=0.0
      a32=0.0
      a33=0.0
      b3=0.0
      dsab=sqrt((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
      if(abs(xa-xb).gt.local_eps*dsab) then
         a21=-(ya-yb)
         a22=  xa-xb
         a23=0.0
         b2=-xb*(ya-yb)+yb*(xa-xb)
         a31=-(za-zb)
         a32=0.0
         a33=  xa-xb
         b3=-xb*(za-zb)+zb*(xa-xb)
      elseif(abs(ya-yb).gt.local_eps*dsab) then
         a21=-(ya-yb)
         a22=  xa-xb
         a23=0.0
         b2=-xb*(ya-yb)+yb*(xa-xb)
         a31=0.0
         a32=-(za-zb)
         a33=  ya-yb
         b3=-yb*(za-zb)+zb*(ya-yb)
      elseif(abs(za-zb).gt.local_eps*dsab) then
         a21=-(za-zb)
         a22=0.0
         a23=  xa-xb
         b2=-xb*(za-zb)+zb*(xa-xb)
         a31=0.0
         a32=-(za-zb)
         a33=  ya-yb
         b3=-yb*(za-zb)+zb*(ya-yb)
      endif
!
!     Normalize the bottom two rows of a and apply the corresponding 
!     transformations to b
      a2norm=sqrt(a21**2+a22**2+a23**2)
      a3norm=sqrt(a31**2+a32**2+a33**2)
!     Row 2
      a21=a21/a2norm
      a22=a22/a2norm
      a23=a23/a2norm
      b2 =b2/a2norm
!     Row 3
      a31=a31/a3norm
      a32=a32/a3norm
      a33=a33/a3norm
      b3 =b3/a3norm
!
!     Make sure the matrix isn't singular, if it is, figure out how far away
!     from the plane the parallel line is. If the line is within releps
!     it's close enough to be considered in the plane
!
!     Calculate the determinant of the line-plane matrix
!
      det = a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+  & 
           a13*(a21*a32-a22*a31)
      if(abs(det).le.local_eps) then
!
!     Find the unit normal to the plane
         lennorm = sqrt(cx**2+cy**2+cz**2)
         xunit = cx/lennorm
         yunit = cy/lennorm
         zunit = cz/lennorm
!     We think we're parallel, make sure...
!     Projection of the line between a point on the plane and a
!     point on the line onto the normal gives us the distance
         dist1 = xunit*(xa-x1)+yunit*(ya-y1)+zunit*(za-z1)
         dist2 = xunit*(xb-x1)+yunit*(yb-y1)+zunit*(zb-z1)
         if(((abs(dist1).le.releps).OR.(abs(dist2).le.releps)) & 
             .AND.(abs(dist1-dist2).le.releps)) then
!     Figure out if either of the two points are within the triangle
!     or not.
            xnodes(1)=x1
            ynodes(1)=y1
            znodes(1)=z1
            xnodes(2)=x2
            ynodes(2)=y2
            znodes(2)=z2
            xnodes(3)=x3
            ynodes(3)=y3
            znodes(3)=z3
!            call inside_element(ifelmtri, xnodes, ynodes, znodes,
!     $           xa, ya, za, iflag)
            call inside_tri( x1, y1, z1, x2, y2, z2, x3, y3, z3, & 
                xa, ya, za, iflag)  
            if(iflag.ge.0) then
               x = xa
               y = ya
               z = za
               iflag = 1
               goto 9999
            endif
!            call inside_element(ifelmtri, xnodes, ynodes, znodes,
!     $           xb, yb, zb, iflag)
            call inside_tri( x1, y1, z1, x2, y2, z2, x3, y3, z3, & 
                xb, yb, zb, iflag)  
            if(iflag.ge.0) then
               x = xb
               y = yb
               z = zb
               iflag = 1
               goto 9999
            endif
!           Ok, the points are not inside the triangle, let's see if the
!           line intersects any of the sides of the triangle
!            do i = 0,2
!               call lineseg_lineseg(xnodes(mod(i,3)+1),
!     &                              ynodes(mod(i,3)+1),
!     &                              znodes(mod(i,3)+1),
!     &                              xnodes(mod(i+1,3)+1),
!     &                              ynodes(mod(i+1,3)+1),
!     &                              znodes(mod(i+1,3)+1),
!     &                              xa,ya,za,
!     &                              xb,yb,zb,
!     &                              iflag)
!               if(iflag.ge.0) then
!                  iflag = 2
!                  goto 9999
!               endif
!            enddo
         elseif(abs(dist1-dist2).lt.releps) then
            iflag = -1
         else
!           Alright, who's making the pathological case?
!           Try and do the best you can with it & warn the user...
            write(warnmess, '(a)')  & 
           'Warning! You may have a pathological case in your problem.'  
            iflag = -1
         endif
         goto 9999
      endif
      iflag=0
!
!     The matrix is not singular, now solve for x, y, z (i.e.,
!     where the plane intersects the line)
!
!     Solve for x (replace a{1,2,3}1 with b{1,2,3})
      x = (b1*(a22*a33-a23*a32)-a12*(b2*a33-a23*b3)+ & 
     &      a13*(b2*a32-a22*b3))/det
!     Solve for y (replace a{1,2,3}2 with b{1,2,3})
      y = (a11*(b2*a33-a23*b3)-b1*(a21*a33-a23*a31)+ & 
     &      a13*(a21*b3-b2*a31))/det
!     Solve for z (replace a{1,2,3}3 with b{1,2,3})
      z = (a11*(a22*b3-b2*a32)-a12*(a21*b3-b2*a31)+  & 
     &      b1*(a21*a32-a22*a31))/det
      if(iflag.eq.0) then
         dsa=sqrt((x-xa)**2+(y-ya)**2+(z-za)**2)
         dsb=sqrt((x-xb)**2+(y-yb)**2+(z-zb)**2)
         dsum=dsa+dsb
         ax123=  (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
         ay123=-((x2-x1)*(z3-z1)-(x3-x1)*(z2-z1))
         az123=  (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
         area123=0.5*sqrt(ax123**2+ay123**2+az123**2)
         ax12=  (y2-y1)*(z-z1)-(y-y1)*(z2-z1)
         ay12=-((x2-x1)*(z-z1)-(x-x1)*(z2-z1))
         az12=  (x2-x1)*(y-y1)-(x-x1)*(y2-y1)
         area12=0.5*sqrt(ax12**2+ay12**2+az12**2)
         ax23=  (y3-y2)*(z-z2)-(y-y2)*(z3-z2)
         ay23=-((x3-x2)*(z-z2)-(x-x2)*(z3-z2))
         az23=  (x3-x2)*(y-y2)-(y3-y2)*(x-x2)
         area23=0.5*sqrt(ax23**2+ay23**2+az23**2)
         ax31=  (y1-y3)*(z-z3)-(y-y3)*(z1-z3)
         ay31=-((x1-x3)*(z-z3)-(x-x3)*(z1-z3))
         az31=  (x1-x3)*(y-y3)-(x-x3)*(y1-y3)
         area31=0.5*sqrt(ax31**2+ay31**2+az31**2)
         xsum=area12+area23+area31
         if(abs(dsab-dsum).le.local_eps*dsab .and.  & 
           abs(area123-xsum).le.local_eps*area123) then
            iflag=0
         else
            iflag=-1
         endif
      endif
      goto 9999
 9999 continue

      pt=(/x,y,z/) 

      d1=sqrt(dot_product(pt-pta,pt-pta) ) 
      d2=sqrt(dot_product(ptb-pta,ptb-pta) ) 

      pt = pta + (ptb-pta) * ( d1/d2 +1.0d-6) 

      return
      end subroutine lineseg_tri
 
 
end module tracktools  
