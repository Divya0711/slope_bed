program genPoints
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,l2
  integer(kind=4)::npt,npm,nx,ny,nz,nth, nDep, nLeftLayer
  integer(kind=4),allocatable::nodid(:), nodEgdeTy(:,:)
  real(kind=8),allocatable::corx(:),cory(:),corz(:)
  real(kind=8),allocatable::SN(:,:),SM(:,:),SS(:,:)
  integer(kind=4)::fsbeg,fsend,wmNod1,wmNod2,nCirLay
  real(kind=8)::domx0,domx1,domy0,domy1,domz0,domz1,domx_f
  real(kind=8)::dx,dy,dz,drMin,pi,tt,dth,slope
  real(kind=8)::cylX,cylY,cylR
  real(kind=8)::tmpr1,tmpr2,tmpr3,tmpr4
  real(kind=8),allocatable::dep(:,:)
  real ab,bb,bb1,s,alpha,Sx,co,co1,cd
	integer(kind=4):: i1,i3,nxx
  
   DOUBLEPRECISION,ALLOCATABLE:: XX(:),HH(:)



  !! nodEgdeTy(a,b)
  !  Edge types defined assuming a rectilinear domain
  !  1,0 = bottom or near edge
  !  2,0 = top or far edge
  !  0,1 = side edge

  allocate(dep(4,2))
  i=1
  dep(i,:) = (/  0.00d0,  0.40d0/); i=i+1;
  dep(i,:) = (/  7.50d0,  0.40d0/); i=i+1;
  dep(i,:) = (/ 13.50d0,  0.10d0/); i=i+1;
  dep(i,:) = (/ 30.00d0,  0.10d0/); i=i+1;
  nDep = i-1

  do k = 1,nDep
    write(*,'(2F15.6)')dep(k,1), dep(k,2)    
  enddo
  write(*,*)

  pi=4d0*atan(1d0)

  domx0=0d0
  domx1=5.0001d0
  domx_f=3.0001d0
  domy0=0d0
  domy1=0.1500001d0
  domz0=0d0
  domz1=0.15000001d0
  dx=0.15d0/20d0
  dy=0.15d0/20d0
  dz=0.15d0/20d0
  drMin=min(dx,dy,dz)/3d0    
  nLeftLayer = 6

  cylX = domx0 + (domx1-domx0)/2d0
  cylY = domy0 + (domy1-domy0)/2d0
  cylR = dx

  nx=floor((domx1-domx0)/dx)
  ny=floor((domy1-domy0)/dy)
  nz=floor((domz1-domz0)/dy)  

  write(*,'(" nx = ",I10)')nx
  write(*,'(" ny = ",I10)')ny
  write(*,'(" nz = ",I10)')nz
  npm=floor((nx+1)*(ny+1)*(nz+1)*1.50)
  write(*,'(" Maxnodes = ",I10)')npm
  npt=0  
  allocate(corx(npm),cory(npm),corz(npm),nodid(-7:npm))
  allocate(SN(npm,3), SM(npm,3), SS(npm,3), nodEgdeTy(npm,2))
  SN=0d0
  SM=0d0
  SS=0d0
  nodEgdeTy=0

  write(*,'(" Check XYZ coord limits and nx,ny,nx limits")')
  write(*,'(" XMax : ",2F15.6)')domx1,(domx0+nx*dx)
  write(*,'(" YMax : ",2F15.6)')domy1,(domy0+ny*dy)
  write(*,'(" ZMax : ",2F15.6)')domz1,(domz0+nz*dz)
  
  !! Fluid particles for slope
    alpha =0.0001
    slope = TAND(5.000d0)
    bb1 = domz1/slope ! getting the base lenght for water depth height
    print*,'slope val',slope,bb1
    pause
    allocate (XX(0:nx+1),HH(0:nx+1))
 
    ab = domx_f
    bb = domx_f+bb1
    nxx = floor(bb/dx) 
    do i =0,nx+1
    xx(i)=i*dx
    end do

    i1 = 0
    i3 = 0

    do j =0,nx+1
    if(xx(j)>=bb) then
    hh(j)=domz1
    elseif ((xx(j)>ab).and.(xx(j)<=bb))then
    hh(j)= domz1-((xx(j)-ab)*((domz1-domz0)/(bb-ab)))
    i3=i3+1
    else
    hh(j)=domz1
    i1=i1+1

    end if

    end do

    deallocate(xx)

     npt =0 
      print*,'nz val',nz

     do k=1,nz
     do j =1,(ny-1)
     do i = nLeftLayer,nxx-1  !nx-1
       tmpr1=domx0+i*dx
        tmpr2=domy0+j*dy  

        npt = npt+1
        s = (alpha*(hh(i))*(nz-k))/(nz-1)
        co = 1-exp(s)
        Sx = alpha*(hh(i))
        co1 = 1-exp(Sx)
        cd =(hh(i))*(co/co1)

         corx(npt)=domx0+i*dx
         cory(npt)=domy0+j*dy
         corz(npt)=domz1-cd
         nodid(npt)=0
     end do
     end do
     end do
      print*,'valu of k after loop',k,nz-1,j,ny-1,i,nxx-1

      print*,'total fluid part',npt
  open(15, file='XY_slp_fluid3d.DAT')
  do i=1,npt
    write(15,'(3F15.6)')corx(i),cory(i),corz(i)
  enddo
    
    print*,'chk fluid'
    pause
    stop
 ! do k=1,nz-1
 !   do j=1,ny-1
  !    do i=-nLeftLayer,nx-1
  !      tmpr1=domx0+i*dx
  !      tmpr2=domy0+j*dy                

   !     npt=npt+1        
   !     corx(npt)=domx0+i*dx
   !     cory(npt)=domy0+j*dy
   !     corz(npt)=domz0+k*dz
    !    nodid(npt)=0      
  !  enddo
  !  enddo
  !enddo  
  nodid(-1)=npt

  !! FS nodes
  fsbeg=npt+1
  do j=1,ny-1
    do i=-nLeftLayer,nx-1
      tmpr1=domx0+i*dx
      tmpr2=domy0+j*dy              

      npt=npt+1
      corx(npt)=domx0+i*dx
      cory(npt)=domy0+j*dy
      corz(npt)=domz1
      nodid(npt)=4
    enddo
  enddo  
  fsend=npt
  nodid(-2)=npt

  !! Wavemaker nodes
  ! do k=0,nz
  !   do j=0,ny
  !     npt=npt+1
  !     corx(npt)=domx0
  !     cory(npt)=domy0+j*dy
  !     corz(npt)=domz0+k*dz
  !     nodid(npt)=8

  !     SN(npt,1) = -1
  !     SM(npt,2) = -1
  !     SS(npt,3) =  1

  !     if(k.eq.0) nodEgdeTy(npt,1) = 1 !bottom
  !     if(k.eq.nz) nodEgdeTy(npt,1) = 2 !top

  !     if((j.eq.0).or.(j.eq.ny)) nodEgdeTy(npt,2)=1 !side

  !     if((k.eq.0).and.(j.eq.floor(ny/2d0))) wmNod1=npt
  !     if((k.eq.nz).and.(j.eq.floor(ny/2d0))) wmNod2=npt
  !   enddo
  ! enddo
  nodid(-3)=npt
  ! write(*,*)corx(wmNod1),cory(wmNod1),corz(wmNod1)
  ! write(*,*)corx(wmNod2),cory(wmNod2),corz(wmNod2)

  !! Bottom particles  
  do j=1,ny-1
    do i=-nLeftLayer,nx-1
      tmpr1=domx0+i*dx
      tmpr2=domy0+j*dy              

      npt=npt+1        
      corx(npt)=domx0+i*dx
      cory(npt)=domy0+j*dy
      corz(npt)=domz0
      nodid(npt)=2

      SN(npt,3) = -1
      SM(npt,2) =  1
      SS(npt,1) =  1
    enddo
  enddo  
  nodid(-4)=npt

  !! Opposite to WM nodes
  do k=0,nz
    do j=0,ny
      npt=npt+1
      corx(npt)=domx1
      cory(npt)=domy0+j*dy
      corz(npt)=domz0+k*dz
      nodid(npt)=3    

      SN(npt,1) =  1
      SM(npt,2) =  1
      SS(npt,3) =  1

      if(k.eq.0) nodEgdeTy(npt,1) = 1 !bottom
      if(k.eq.nz) nodEgdeTy(npt,1) = 2 !top

      if((j.eq.0).or.(j.eq.ny)) nodEgdeTy(npt,2)=1 !side

    enddo
  enddo
  nodid(-5)=npt

  !! Near Side wall
  do k=0,nz
    do i=-nLeftLayer,nx-1
      npt=npt+1
      corx(npt)=domx0+i*dx
      cory(npt)=domy0
      corz(npt)=domz0+k*dz
      nodid(npt)=1

      SN(npt,2) = -1
      SM(npt,1) =  1
      SS(npt,3) =  1

      if(k.eq.0) nodEgdeTy(npt,1) = 1 !bottom
      if(k.eq.nz) nodEgdeTy(npt,1) = 2 !top
      
    enddo
  enddo
  nodid(-6)=npt

  !! Far Side wall
  do k=0,nz
    do i=-nLeftLayer,nx-1
      npt=npt+1
      corx(npt)=domx0+i*dx
      cory(npt)=domy1
      corz(npt)=domz0+k*dz
      nodid(npt)=7

      SN(npt,2) =  1
      SM(npt,1) = -1
      SS(npt,3) =  1

      if(k.eq.0) nodEgdeTy(npt,1) = 1 !bottom
      if(k.eq.nz) nodEgdeTy(npt,1) = 2 !top
      
    enddo
  enddo
  nodid(-7)=npt

  !! Cylinder wall
  npt=npt+1
  corx(npt)=cylX+cylR
  cory(npt)=cylY
  corz(npt)=2d0*domz1
  nodid(npt)=9
  SN(npt,1) = -1
  SM(npt,2) = -1
  SS(npt,3) =  1

  npt=npt+1
  corx(npt)=cylX-cylR
  cory(npt)=cylY
  corz(npt)=2d0*domz1
  nodid(npt)=9
  SN(npt,1) =  1
  SM(npt,2) =  1
  SS(npt,3) =  1

  npt=npt+1
  corx(npt)=cylX
  cory(npt)=cylY+cylR
  corz(npt)=2d0*domz1
  nodid(npt)=9
  SN(npt,2) = -1
  SM(npt,1) =  1
  SS(npt,3) =  1

  npt=npt+1
  corx(npt)=cylX
  cory(npt)=cylY-cylR
  corz(npt)=2d0*domz1
  nodid(npt)=9
  SN(npt,2) =  1
  SM(npt,1) = -1
  SS(npt,3) =  1
  
  nodid(-1)=npt

  !! Writing the mesh file
  open(11, file='mesh.dat')
  write(11,'("NODEID(-1:-7)")')  
  write(11,'(7I15)')(nodid(j),j=-1,-7,-1)
  write(11,'("Free Surface Indices")')  
  write(11,'(2I15)')fsbeg,fsend
  write(11,'("Wavemaker reference nodes & their coors")')  
  write(11,'(I15,3F15.6)')wmNod1,corx(wmNod1),&
    cory(wmNod1),corz(wmNod1)
  write(11,'(I15,3F15.6)')wmNod2,corx(wmNod2),&
    cory(wmNod2),corz(wmNod2)
  write(11,'("Node Coordinates")')      
  do i=1,npt
    write(11,'(12F15.6, 2I4)')corx(i), cory(i), corz(i), &
      SN(i,1:3), SM(i,1:3), SS(i,1:3), nodEgdeTy(i,1:2)
  enddo
  close(11)

  !! Writing output file for visualisation in Paraview
  open(12, file='XY_LOCATION.DAT')
  write(12,'(I18)')npt
  write(12,*)'TIME=',0D0,0
  do i = 1, npt
    write(12,'(8F15.6,1I4)')corx(i),cory(i),corz(i),&
      0D0,SN(i,1:3),0d0,nodid(i)
  enddo
  write(12,*)'***'
  close(12)

end program




subroutine interpDep(nDep, dep, x, d)
implicit none
  
  integer(kind=4),intent(in):: nDep
  real(kind=8),intent(in):: dep(nDep,2), x
  real(kind=8),intent(out):: d

  integer(kind=4):: i
  real(kind=8):: sl

  if((x.lt.dep(1,1)).or.(x.gt.dep(nDep,1)))then
    write(*,*)"Error in dep range for x", x
    stop
  endif

  do i = 1, nDep
    if(dep(i,1).gt.x)exit
  enddo

  sl = (dep(i,2)-dep(i-1,2))/(dep(i,1)-dep(i-1,1))
  d = dep(i-1,2) + sl*(x - dep(i-1,1))

end subroutine interpDep
