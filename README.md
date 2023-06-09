# slope_bed
generation of slope bed in 3d imlpgr

the mesh gen file from shagun is modified to incorporate slope bed with left buffer layer to capture fnpt variables

- the fluid particle generation is done for slope bed. The do loop 
   { alpha =0.0001
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
    }
    
    is added to  get the depth value for different X 
    
    - here domz1 value is reduced to one layer below the freesurface layer so that the free surface nodes are stored separately for future purpose
    - domz1=(0.15000001d0-dz) is included for this reason so that the fluid particles are accomodated from one layer below fsf for sloping bed
