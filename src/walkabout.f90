program walkabout  

!******************************************************************
! Walkabout 
! Random walk particle tracking on unstructured finite volume grids 
! Developed by Scott Painter 
! Version 1.2
! March 2011 
! 
! This material was prepared by Los Alamos National Security, LLC (LANS) 
! under Contract DE-AC52-06NA25396 with the U.S. Department of Energy 
! (DOE). All rights in the material are reserved by DOE on behalf of the 
! Government and LANS pursuant to the contract. You are authorized to use 
! the material for Government purposes but it is not to be released or 
! distributed to the public. NEITHER THE UNITED STATES NOR THE UNITED 
! STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC, 
! NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
! ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, 
! COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR 
! PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE 
! PRIVATELY OWNED RIGHTS.
!******************************************************************  

!Parallelized Walkabout using OpenMP
!Implemented by Soumi Manna
!September, 2015

use precision 
use iotools 
use mesh_class 
use cell_class 
use flowfield_class 
use tracktools 
use omp_lib

implicit none 

character(len=80) :: ptitle ! title of simulation. 1st line in control.dat 
character(len=70) :: filesfile, outfile , sptr2file 
character(len=70) :: finposfile 
integer :: nn,m,n , ncells,melem , istep ,  nstep ,icol , maxsteps 
integer :: firstelem 
integer :: isolution  , itime , nsolutions, numargs 
integer :: i,j,k,l 
integer :: AllocateStatus
integer :: ireject, nreject=50  ! keep track of rejected dispersive steps 

logical :: rejected = .false., noflow = .false. 
integer :: iprog 
real(kind=dp) :: nextprog  
real(kind=dp), dimension(7), parameter :: prog = (/ 0.05, 0.10, 0.20, 0.50, 0.7, 0.9, 1.01 /)  

integer :: oldcell, newcell 
integer :: ipart, npart 
real(kind=dp), dimension(:,:), allocatable :: starts 
real ( kind = 8 ) :: time_begin
real(kind = 8) :: time_end
real ( kind = 8 ) :: time_total
logical :: plottraj ! flag to plot trajectories  
integer :: ptfreq ! frequency of trajectory output  
integer ::  numout ! number trajectory outputs each particle 

type(cell), pointer :: ptr2cell ,ptr2nbr 
type(element), pointer :: ptr2elem 
type(mesh), target  :: crntmesh 
type(mesh), pointer :: ptr2mesh 

real(kind=dp) ::  time , dxtar, dttar,  maxstretch 
real(kind=8), dimension(:,:), allocatable :: posrecord
integer, dimension(:), allocatable :: numout_part
real(kind=dp) ::  dt, dt0, dt1, dt2, dtmax
type(flowfield), dimension(:), allocatable,target :: fehmresults 
type(flowfield), pointer :: crntfehm  

real(kind=dp),dimension(3) :: vel ,  pt , dx , workpt , delv , crosspt 
real(kind=dp),dimension(3) :: vel1 , disp, trial  , drift 
real(kind=dp), dimension(3,3) :: dlocal 

real(kind=dp), dimension(3) :: vert1, vert2, vert3, vert4 
real(kind=dp), dimension(4,3) :: vertvels 
real(kind=dp) :: z1, z2, z3 , absvel 
real(kind=dp) :: partdisp, totdisp , ptime_next
real :: t0, t1, t2, t3, t4

real(kind=dp), parameter :: ascale=sqrt(3.0d0) 
integer :: nextelem  , trialelem 
integer,dimension(12) :: seed
integer :: seedsize 
integer :: sunit,tunit,tmpunit, fpunit,chunk,ichunk,nchunk 
logical :: enableNID , bflag 
bflag = .false.

! print the version number 
call printversion() 

! initialize cpu timers 
call cpu_time(t0) 

! get names file from command line 
numargs = iargc() 
if(numargs .eq. 0) then   
 filesfile="walkabout.files"  ! default 
else 
 call getarg(1, filesfile) 
end if 


! the following will load various control parameters  
! also sets output file and some output options 
call getcontrols() 

! new get starting locations for all particles (npart, and start array)  
call getstarts() 

!obtain processor and thread number

write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

! initialize random seeds 
!  if seed is not already set as input then initialize from clock   
! initialize random seeds
!  if seed is not already set as input then initialize from clock
if( seed(1) .eq. 0) then
  call date_and_time(VALUES=seed)          ! get the time of day
  call random_seed(PUT=seed(1:12:1) )      ! initialize the seed based on time and date
else
  call random_seed(put=seed)
end if

! load the mesh 
crntmesh = mesh_( filesfile )
ptr2mesh => crntmesh 

call cpu_time(t1) 

print *, 'successfully read ' , crntmesh%nelems, ' elements' 
print *, 'successfully read ', crntmesh%ncells, ' cells' 
print *, ' time= ', t1-t0 

  ! will need to allocate array of solution type to hold time-dependent results 
  nsolutions=1 
  allocate(fehmresults(nsolutions)) 
  fehmresults(1)=flowfield_(filesfile,ptr2mesh)  
  ! could put in loop below then release the flux array when no longer needed 


! for each solution from FEHM 

  do isolution=1,size(fehmresults) 
    crntfehm => fehmresults(isolution) 
    call reconstruct_flowfield( crntfehm) 
  end do ! velocity fields for all times now loaded 

  call cpu_time(t2) 
  print *, 'finished velocity reconstruction' 
  print *, ' time= ', t2-t1 

! ************************
! tracking starts here 

  sunit=findunitnum() 
  open(unit=sunit, file=sptr2file) 
  call printheader( sunit, ptitle ) 
  write(sunit,*) npart 
  write(sunit,*) ' Part_no    time_days     cell_leaving' 

  fpunit=findunitnum() 
  open(unit=fpunit, file=finposfile) 
  call printheader( fpunit, ptitle ) 
  write(fpunit,*) npart 
  
    
  iprog=1 
  nextprog=prog(iprog) 
  ptr2mesh => crntmesh 
  pt=starts(1,:) 
  !firstelem=spiralsearch(ptr2mesh, 1, pt)
  firstelem=findelement(ptr2mesh, pt)
  !For border particle release debugging
  print *, "we found the first element at: ", firstelem
  
  if (ptfreq .eq. 0) then
  	  nonzero_freq = 1
  else
  	  nonzero_freq = ptfreq
  endif 
  
  !write(*,*) 1+maxsteps/nonzero_freq

tunit=findunitnum()
 if(plottraj) then
    open(unit=tunit, file=outfile)
   call printheader(tunit, ptitle)
    write(tunit,*) npart 
endif

  !Measure total time before the OpenMP loop parallelization starts
  time_total=0
  chunk = 10000
  if(npart .lt. chunk) then
    chunk = npart
  endif
     
  ! Number of particles are divided by chunks and provide better parallel performance
  nchunk = npart / chunk
  
  !Temporary array for reading and writing particles tracks
  ALLOCATE(posrecord(2+(maxsteps/nonzero_freq),chunk*4))
  allocate(numout_part(chunk))
  
   do ichunk=1, nchunk

   !Add OpenMP loop parallelization with private and shared data

   !$omp parallel default(shared) private(time,dt,pt,numout, bflag,melem,crntfehm,istep,ptr2elem, oldcell, itime,  enableNID, vel,dt1,trial,drift,nextelem,workpt, vel1,dx,rejected,ireject,trialelem,i,j,k,l, vert1,vert2,vert3,vert4,icol,dlocal, disp, newcell, vertvels)


    !Calculate time at the begining of the openmp directive
    time_begin = omp_get_wtime()
     
    !$omp do schedule(dynamic)
    
    do ipart=1,chunk 
    
    time=0.0 

    pt=starts(chunk*(ichunk-1)+ipart,:) 
     
    
    if( plottraj) then
       numout_part(ipart)=1
       posrecord(numout_part(ipart),(ipart-1)*4+1) = time
       posrecord(numout_part(ipart),(ipart-1)*4+2) = pt(1)
       posrecord(numout_part(ipart),(ipart-1)*4+3) = pt(2)
       posrecord(numout_part(ipart),(ipart-1)*4+4) = pt(3)
       
       numout_part(ipart)=numout_part(ipart)+1
    endif
    
    dt = dt0  
    melem=spiralsearch(ptr2mesh, firstelem, pt) 
    ptr2elem => crntmesh%elems(melem)
    
    oldcell = findcell( ptr2mesh, ptr2elem, pt)  

      !particleloop: 
              
            do istep=1,maxsteps 
            ! find the flow field active at current time 
            do itime=size(fehmresults) , 1, -1 
               if( time .ge. fehmresults(itime)%time ) exit  
            end do  
            crntfehm => fehmresults(itime)   !  flow field applicable at given time  

            ! effective velocity including NID 
            ptr2elem => crntmesh%elems(melem) 
            enableNID=.true. 
            vel=veffective(ptr2elem, ptr2mesh, crntfehm, pt, enableNID) 

            ! time step control 
            dt1= dtlimit( ptr2elem, ptr2mesh, crntfehm, pt, vel, dxtar, dttar)
            dt = min( dt1, dtmax, dt*maxstretch) 
            
            ! predictor 
            1001 drift=dt*vel 
             trial=pt+drift 

            ! corrector step 
            !  disable noise-induced drift if corrector step crosses a boundary 
            !  advective component will extrapolate in that situation 
            !  scheme prevents drift across a no-flow boundary 
            nextelem = melem 
            workpt=pt 
            call intersect_face(ptr2mesh, nextelem, workpt, trial, bflag, noflow) 
            if( bflag ) then  
              trial=workpt 
              enableNID=.false. 
            endif 

            ptr2elem => crntmesh%elems(nextelem)  
            vel1=veffective(ptr2elem, ptr2mesh, crntfehm, trial,enableNID) 

            ! time step control v1.2 
            !  reject time step and try again if needed 
            dt1= dtlimit( ptr2elem, ptr2mesh, crntfehm, trial, vel1, dxtar, dttar)
            if( dt1 .lt. dt/2.0d0) then 
              dt=dt/2.0d0 
              goto 1001 
            end if 

            dx=vel1*dt 
            trial=pt + dx  ! could weight between predictor and corrector? 

            nextelem=melem 
            workpt=pt 
            call intersect_face(ptr2mesh, nextelem, workpt, trial, bflag, noflow) 
            ! cut time step near no-flow boundary
            if( noflow .and. .not.bflag) then 
              print *, ' cannot happen. error: '
              stop 
            end if 
            if( noflow ) then 
              dt=dt/2.0d0 
              goto 1001 
            end if 
            if(bflag) exit 

            trialelem=nextelem  !element for dispersive step 

            ! advection/drift step finished now disperse 
            ! only complication here is that particle is not allowed to disperse across
            !   a zero concentration boundary 

            rejected = .true. 
            ireject=0  
            do while( rejected .and. ireject .lt. nreject) 

                ! random numbers 
                ! ascale makes random variable have unit variance  
                call random_number(z1) 
                call random_number(z2) 
                call random_number(z3) 
                z1=(2.0d0*z1-1.0d0)*ascale 
                z2=(2.0d0*z2-1.0d0)*ascale 
                z3=(2.0d0*z3-1.0d0)*ascale 

                ! interpolate dispersive displacement tensor, one column at a time 
                ! melem and pt are prior to advective step, consistent with Itoh's SDE formulation 
                ptr2elem => crntmesh%elems(melem) 
                i=ptr2elem%i
                j=ptr2elem%j
                k=ptr2elem%k
                l=ptr2elem%l
                vert1 = crntmesh%cells( i )%posi
                vert2 = crntmesh%cells( j )%posi
                vert3 = crntmesh%cells( k )%posi
                vert4 = crntmesh%cells( l )%posi
               
                do icol=1,3 
                    vertvels(1,:) = crntfehm%disp(i,:,icol)
                    vertvels(2,:) = crntfehm%disp(j,:,icol)
                    vertvels(3,:) = crntfehm%disp(k,:,icol)
                    vertvels(4,:) = crntfehm%disp(l,:,icol)
                    dlocal(:,icol) =  interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt,0)
                end do
                
 
                ! dispersive displacement 
                disp(1) = (z1*dlocal(1,1) + z2*dlocal(1,2) + z3*dlocal(1,3) ) * sqrt(dt)  
                disp(2) = (z1*dlocal(2,1) + z2*dlocal(2,2) + z3*dlocal(2,3) ) * sqrt(dt)  
                disp(3) = (z1*dlocal(3,1) + z2*dlocal(3,2) + z3*dlocal(3,3) ) * sqrt(dt)  

                workpt=trial 
                nextelem = trialelem 
                call intersect_face(ptr2mesh, nextelem, workpt, trial+disp, bflag, rejected) 
    
                if( rejected ) ireject=ireject+1 
 
            end do  

            pt = trial + disp
            
            time=time+dt
             
            
             if(plottraj .and.  mod(istep,ptfreq) .eq. 0 ) then

              
              posrecord(numout_part(ipart),(ipart-1)*4+1) = time
              posrecord(numout_part(ipart),(ipart-1)*4+2) = pt(1)
              posrecord(numout_part(ipart),(ipart-1)*4+3) = pt(2)
              posrecord(numout_part(ipart),(ipart-1)*4+4) = pt(3)
              numout_part(ipart)=numout_part(ipart)+1
               
            end if 
            

            if( bflag) exit  

            melem=nextelem 
    
            ptr2elem => crntmesh%elems(melem) 
            newcell = findcell( ptr2mesh, ptr2elem, pt)  
            if( oldcell .ne. newcell) then 
               write(sunit,'(i9,2x,g15.8,i11)'), ipart, time, oldcell
            end if 
            oldcell=newcell 
 
        end do !particleloop
        
        if( npart .gt. 20) then 
            if(  ipart .gt. nextprog*npart) then 
                write(*,'(i2,a10)') , int( nextprog*100) , '% complete'
                if (iprog < 7) then
                	iprog=iprog+1      
                endif
                nextprog = prog(iprog) 
            end if 
        end if 
        ! write final position if plottraj selected 
        if(plottraj) then
              
              
              posrecord(numout_part(ipart),(ipart-1)*4+1) = time
              posrecord(numout_part(ipart),(ipart-1)*4+2) = pt(1)
              posrecord(numout_part(ipart),(ipart-1)*4+3) = pt(2)
              posrecord(numout_part(ipart),(ipart-1)*4+4) = pt(3)
              
        end if 

        write(sunit,'(i9,2x,g15.8,i11)'), ipart, time, oldcell

        !write final position to finalposition file 
        write(fpunit,'(4g15.8)' ) time,pt 
        
    

end do ! end of all particles

!!$omp end do nowait 
!$omp end parallel

!End of OpenMP loop parallelization

!Calculate end time after the OpenMP directive
time_end = omp_get_wtime()

!Calculate total time with start and end time of the OpenMP directive
time_total=time_total+ time_end-time_begin

! Write out posrecord (store the particle tracks)
 do ipart=1, chunk
    nstep=istep
    write(tunit,*) numout_part(ipart)
    do istep=1,numout_part(ipart) 
       write(tunit,'(4g15.8)' ) posrecord(istep,(ipart-1)*4+1),posrecord(istep,(ipart-1)*4+2:(ipart-1)*4+4)
    end do
   
       !write(*,*) chunk*(ichunk-1)+ipart, 'finished'
       
 enddo

end do !chunk loop
close(tunit)

!Print the total time
print *, 'Time of operation was ', time_total, ' seconds'
 
    close(sunit) 
    close(fpunit) 


    call cpu_time(t3) 
    print *, 'cpu time for particle tracking: ', t3-t2 
    print *, 'cpu time for velocity reconstruction ', t2-t1
    print *, 'cpu time for reading ', t1-t0 
    print *, 'total cpu time ', t3-t0 
    print *, ' normal end walkabout' 
!endif 
stop 

contains 

!*********************************************************************************
subroutine getcontrols() 
!  module procedure
!  loads control parameters 

 character(len=80) :: buffer 
 character(len=70) :: controlfile  
 character(len=30) :: key 
 integer :: funit 

 funit=findunitnum() 
 open(unit=funit,file=filesfile,action='read')

 ! file for trajectories 
 outfile="walkabout.out"
 do
  read(unit=funit,fmt='(a)', end=999) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8) .eq. "trajout:") exit
 end do
 outfile = buffer(9:)
 print *, 'trajout file: ', outfile 
999 continue 

 ! sptr2 file  
 rewind(funit) 
 sptr2file="walkabout.sptr2"
 do
  read(unit=funit,fmt='(a)', end=998) buffer
  buffer=adjustl(buffer)
  if( buffer(1:6) .eq. "sptr2:") exit
 end do
 sptr2file = buffer(7:)
998 continue 
 print *, 'sptr2 format for plumecalc: ', sptr2file  

 ! control file 
 rewind(funit) 
 do
  read(unit=funit,fmt='(a)', end=997) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8) .eq. "control:") exit
 end do
 controlfile = buffer(9:)
 print *, 'control file: ', controlfile  

 ! file for final positions 
 rewind(funit) 
 finposfile="walkabout.finpos" 
 do
  read(unit=funit,fmt='(a)', end=996) buffer
  buffer=adjustl(buffer)
  if( buffer(1:7) .eq. "finpos:") exit
 end do
 finposfile = buffer(8:)
996 continue 
 print *, 'final position file: ', finposfile 

 close(funit) 
 
 ! now open control file and read control parameters 
  
 funit=findunitnum() 
 open(unit=funit, file=controlfile, action='read') 

 ! ptitle 
 read(unit=funit,fmt='(a)') ptitle  

 key='dt0'
 dt0=0.01  ! default 
 do
  read(unit=funit,fmt='(a)', end=99) buffer
  buffer=adjustl(buffer)
  if( buffer(1:3).eq.key) then 
    read(buffer(5:), *) dt0  
    exit 
  end if 
 end do 
 99 continue  

 rewind(funit) 
 key='dxtarget' 
 dxtar=0.1  ! relative to cell size  
 do
  read(unit=funit,fmt='(a)', end=98) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8).eq.key) then 
    read(buffer(9:), *) dxtar  
    exit 
  end if 
 end do 
 98 continue  

 rewind(funit) 
 key='maxstretch' 
 maxstretch=1.2  ! default 
 do
  read(unit=funit,fmt='(a)', end=97) buffer
  buffer=adjustl(buffer)
  if( buffer(1:10).eq.key) then 
    read(buffer(11:), *) maxstretch  
    exit 
  end if 
 end do 
 97 continue  

 rewind(funit) 
 key='dtmax' !days 
 dtmax=100000.  ! default 
 do
  read(unit=funit,fmt='(a)', end=96) buffer
  buffer=adjustl(buffer)
  if( buffer(1:5).eq.key) then 
    read(buffer(6:), *) dtmax  
    exit 
  end if 
 end do 
 96 continue  

 rewind(funit) 
 key='maxsteps' !days 
 maxsteps=100000  ! default 
 do
  read(unit=funit,fmt='(a)', end=95) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8).eq.key) then 
    read(buffer(9:), *) maxsteps  
    exit 
  end if 
 end do 
 95 continue  

 rewind(funit) 
 key='dttarget'
 dttar=0.1  ! default relative to characteristic time for dispersion 
 do
  read(unit=funit,fmt='(a)', end=94) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8).eq.key) then
    read(buffer(9:), *) dttar 
    exit
  end if
 end do
 94 continue

 rewind(funit) 
 key='toutfreq' 
 ptfreq=0  ! default , no traj.out output
 do
  read(unit=funit,fmt='(a)', end=93) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8).eq.key) then
    read(buffer(9:), *) ptfreq 
    exit
  end if
 end do
 93 continue
 plottraj=.false.
 if( ptfreq .gt. 0) plottraj=.true.  ! not trajectory output 

 rewind(funit) 
 key='seed' 
 seed=0  
 do
  read(unit=funit,fmt='(a)', end=92) buffer
  buffer=adjustl(buffer)
  if( buffer(1:4).eq.key) then
    read(buffer(5:), *) seed(1)  
    exit
  end if
 end do
 92 continue
  

 print *, 'particle tracking control parameters' 
 print *, 'dt0: ', dt0 
 print *, 'dxtarget: ', dxtar 
 print *, 'dttarget: ', dttar 
 print *, 'maxstretch: ', maxstretch 
 print *, 'dtmax: ', dtmax 
 print *, 'maxsteps: ', maxsteps 
 print *, 'toutfreq: ', ptfreq 
 print *, 'seed: ', seed(1) 
 print *, ' ' 

 close(funit) 
 return 

 997 print *, 'control filename not specified'
 stop 
end subroutine getcontrols 

!****************************************************************
subroutine getstarts 
! gets starting positions 
! module procedure 

 character(len=80) :: buffer 
 character(len=70) :: controlfile  
 character(len=30) :: key 
 real(kind=dp) :: r1,r2,r3 
 integer :: i, j, k , nx, ny, nz, values(1:8) 
 real(kind=dp) :: x1, x2, y1, y2, z1, z2, dx, dy, dz 
 real(kind=dp) :: xpos, ypos, zpos 
 integer, allocatable, dimension(:) :: seed
 integer :: funit, iseed 

    open(52,file='initial_pos.dat')

 funit=findunitnum() 
 open(unit=funit,file=filesfile,action='read') 
 do
  read(unit=funit,fmt='(a)', end=99) buffer
  buffer=adjustl(buffer)
  if( buffer(1:8) .eq. "control:") exit
 end do
 controlfile = buffer(9:)
 close(funit) 
 
 key='INITIAL' 
 funit=findunitnum() 
 open(unit=funit, file=controlfile, action='read') 
 do
  read(unit=funit,fmt='(a)', end=99) buffer
  buffer=adjustl(buffer)
  buffer=trim(buffer)
  if( buffer(1:7) .eq. 'INITIAL') exit 
 end do 

 read(unit=funit,fmt='(a)') buffer
 buffer=adjustl(buffer)
 select case(buffer(1:6)) 
  case('MANUAL')  
   read(unit=funit, fmt=*) npart  
   allocate( starts(npart,3) ) 
   do ipart=1,npart 
     read(unit=funit, fmt=*) starts(ipart,:) 
   end do 
  case('RANDOM') 
    read(unit=funit, fmt=*) npart  
    allocate( starts(npart,3) ) 
    read(unit=funit, fmt=*) x1, y1, z1 
    read(unit=funit, fmt=*) x2, y2, z2 

!   This following 5 lines are added to create different random numbers 
!   in different runs.  Zhiming Lu (2/3/2015)
    call date_and_time(values=values) 
    CALL RANDOM_SEED(size=iseed)
    allocate(seed(1:iseed))
    seed(:) = values(8)
    call random_seed(put=seed)

    do ipart=1,npart 
      call random_number(r1) 
      call random_number(r2) 
      call random_number(r3) 
      starts(ipart,1)=x1+(x2-x1)*r1 
      starts(ipart,2)=y1+(y2-y1)*r2 
      starts(ipart,3)=z1+(z2-z1)*r3 
    end do 
  case('UNIFOR')  
    read(unit=funit, fmt=*) nx, ny, nz 
    npart=nx*ny*nz 
    allocate( starts(npart,3) ) 
    read(unit=funit, fmt=*) x1, y1, z1 
    read(unit=funit, fmt=*) x2, y2, z2 
    dx=(x2-x1)/float(nx) 
    dy=(y2-y1)/float(ny) 
    dz=(z2-z1)/float(nz) 
    ipart=0
    do i=1,nx 
     xpos = x1 + dx/2 + (i-1)*dx 
     do j=1,ny
       ypos = y1 + dy/2 + (j-1)*dy 
       do k=1,nz 
         zpos = z1 + dz/2 + (k-1)*dz 
         ipart=ipart + 1 
         starts(ipart, 1) = xpos 
         starts(ipart, 2) = ypos 
         starts(ipart, 3) = zpos 
       end do 
     end do 
    end do 
  case DEFAULT 
    print *, ' must specify starting positions' 
    stop 
 end select 

! Z. Lu
!  write initial positions for checking results
   do ipart=1,npart
      write(52,'(i10,3f15.4)') ipart, starts(ipart,1), starts(ipart,2), starts(ipart,3)
   enddo
   close(52)

 close(funit)   
 return 

99 print *, ' no initial positions specified. Stopping' 
  stop 

end subroutine getstarts 

!*****************************************************************
function dtlimit( ptr2elem, crntmesh, crntfehm, pt,vel,dxtar,dttar) result(dt) 

type(element), pointer :: ptr2elem
type(mesh), pointer :: crntmesh
type(flowfield), pointer :: crntfehm
real(kind=dp) :: dt1, dt2 , dt, v , dxtar, dttar  , lbar 
real(kind=dp), dimension(3) :: pt,vel 
real(kind=dp), dimension(3) ::  vert1, vert2, vert3, vert4
real(kind=dp), dimension(4) :: vertvels
real(kind=dp),dimension(:),pointer :: dtdisp  
real(kind=dp),dimension(:),pointer :: vols 
integer :: i,j,k,l

    dtdisp => crntfehm%dtdisp 
    vols => crntmesh%vols  

    i=ptr2elem%i
    j=ptr2elem%j
    k=ptr2elem%k
    l=ptr2elem%l
    vert1 = crntmesh%cells( i )%posi
    vert2 = crntmesh%cells( j )%posi
    vert3 = crntmesh%cells( k )%posi
    vert4 = crntmesh%cells( l )%posi
    vertvels(1) = dtdisp(i) 
    vertvels(2) = dtdisp(j)  
    vertvels(3) = dtdisp(k)  
    vertvels(4) = dtdisp(l)  
    dt1=interpolate_tet1d( vert1, vert2, vert3, vert4, vertvels, pt,0) *dttar 

    vertvels(1) = vols(i) 
    vertvels(2) = vols(j)  
    vertvels(3) = vols(k)  
    vertvels(4) = vols(l)  
    v=interpolate_tet1d( vert1, vert2, vert3, vert4, vertvels, pt,0) 
    lbar = v**(1.0d0/3.0d0)
    dt2 = lbar/sqrt(dot_product(vel,vel))  *dxtar 

    dt=min(dt1,dt2) 
!   add by ZL for checking dt
    if( vertvels(1) < 0.0d0 .or. vertvels(2) < 0.0d0 .or. vertvels(3) < 0.0d0 .or. vertvels(4) < 0.0d0) dt =dt0

    return 

end function dtlimit 
!*******************************************************************
function veffective( ptr2elem, crntmesh, crntfehm, pt, enableNID) result(vel) 

type(element), pointer :: ptr2elem
type(mesh), pointer :: crntmesh 
type(flowfield), pointer :: crntfehm 
logical :: enableNID 
real(kind=dp), dimension(3) :: pt , gradpor, vel1 
real(kind=dp) :: por 
real(kind=dp), dimension(3) :: vel, vert1, vert2, vert3, vert4  
real(kind=dp), dimension(3,3) :: dlocal 
real(kind=dp), dimension(4,3) :: vertvels  
real(kind=dp), dimension(4) :: vals 
real(kind=dp), dimension(:,:), pointer :: vfield  
real(kind=dp), dimension(:,:,:), pointer :: dcoef 
real(kind=dp), dimension(:), pointer :: porosity 
integer :: i,j,k,l 

vfield => crntfehm%vfield 
dcoef => crntfehm%dcoef 
porosity => crntfehm%porosity  

    i=ptr2elem%i
    j=ptr2elem%j
    k=ptr2elem%k
    l=ptr2elem%l
    vert1 = crntmesh%cells( i )%posi
    vert2 = crntmesh%cells( j )%posi
    vert3 = crntmesh%cells( k )%posi
    vert4 = crntmesh%cells( l )%posi
    vertvels(1,:) = vfield(i,:)
    vertvels(2,:) = vfield(j,:)
    vertvels(3,:) = vfield(k,:)
    vertvels(4,:) = vfield(l,:)
    vel=interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 0)

    if(.not.enableNID) return  

    vertvels(1,:) = dcoef(i,:,1)
    vertvels(2,:) = dcoef(j,:,1)
    vertvels(3,:) = dcoef(k,:,1)
    vertvels(4,:) = dcoef(l,:,1)
    vel=vel + interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 1)
    dlocal(:,1) = interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 0)

    vertvels(1,:) = dcoef(i,:,2)
    vertvels(2,:) = dcoef(j,:,2)
    vertvels(3,:) = dcoef(k,:,2)
    vertvels(4,:) = dcoef(l,:,2)
    vel=vel + interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 2)
    dlocal(:,2) = interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 0)

    vertvels(1,:) = dcoef(i,:,3)
    vertvels(2,:) = dcoef(j,:,3)
    vertvels(3,:) = dcoef(k,:,3)
    vertvels(4,:) = dcoef(l,:,3)
    vel=vel + interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 3)
    dlocal(:,3) = interpolate_tet3d( vert1, vert2, vert3, vert4, vertvels, pt, 0)

    if(.not. crntfehm%amaread)then
        ! porosity drift term 
        vals(1) = porosity(i)
        vals(2) = porosity(j)
        vals(3) = porosity(k)
        vals(4) = porosity(l)
        por= interpolate_tet1d( vert1, vert2, vert3, vert4, vals, pt, 0)
        gradpor(1)= interpolate_tet1d( vert1, vert2, vert3, vert4, vals, pt, 1)
        gradpor(2)= interpolate_tet1d( vert1, vert2, vert3, vert4, vals, pt, 2)
        gradpor(3)= interpolate_tet1d( vert1, vert2, vert3, vert4, vals, pt, 3)

        do i=1,3 
         vel1(i) = dot_product( dlocal(i,:) , gradpor ) 
        end do 
        vel1 = vel1/por 

        vel=vel + vel1  
    endif

    return 

end function veffective 

end program walkabout  

