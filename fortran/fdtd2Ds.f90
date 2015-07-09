program fdtd2Ds

  implicit none

  integer :: nx,ny,nx2,ny2,isx,isy,irx,iry,nt,npml,pmlfac,pmlexp
  real    :: vp0, rho, fc, dx, dy, dt, tsour, lsour, t0, a
  real, allocatable :: t(:), fs(:), sfd(:), tau(:)
  real, allocatable :: qx(:,:),qy(:,:),ux(:,:),uy(:,:),px(:,:),py(:,:),vp(:,:),p(:,:)
  !real, allocatable :: diffop(:,:), pmlop(:,:)
  real:: diffop, pmlop
  real, parameter :: pi = acos(-1.)
  character(len=200) :: number
  integer :: i,j,k,l

  ! define velocity model
  nx = 400
  ny = 400
  nx2 = nx+1
  ny2 = ny+1
  vp0 = 2000.

  allocate(vp(1:ny2,1:nx2))
  vp = vp0**2
  rho = 1500.

  ! define source position and receiver position
  isx = 200
  isy = 100
  irx = 200
  iry = 300

  ! source and fd parameters
  nt = 3000
  fc = 1000.
  dx = vp0/(30*fc)
  dy = dx
  dt = dx/vp0*1/(sqrt(2.)*2*pi)

  ! define source function in time
  tsour = 1/fc
  allocate(t(nt))
  t=0.
  t(1:nt) =  (/ ((i-1)*dt, i=1,nt) /) ! implicit loop
  lsour = tsour/dt
  t0 = tsour*1.5
  allocate(tau(1:nt))
  tau = pi*(t-t0)/t0
  a=4
  allocate(fs(1:nt))
  fs=0.
  fs = (1-a*tau*tau)*exp(-2*tau*tau)

  ! define PML properties
  npml = 30
  pmlfac = 50
  pmlexp = 2
  allocate(qx(1:ny2,1:nx2))
  allocate(qy(1:ny2,1:nx2))
  qx = 0.
  qy = 0.
  do i=1,npml
    qx(1:ny2,i) = float(pmlfac*(npml-i)**pmlexp) !left
    qx(1:ny2,nx2-i+1) = float(pmlfac*(npml-i)**pmlexp) !rmght
    qy(i,1:nx2) = float(pmlfac*(npml-i)**pmlexp) !top
    qy(ny2-i+1,1:nx2) = float(pmlfac*(npml-i)**pmlexp) !bottom
  enddo

  !initialize fields
  allocate(px(1:ny2,1:nx2), py(1:ny2,1:nx2)) !pressure field
  allocate(ux(1:ny2,1:nx2), uy(1:ny2,1:nx2)) !particle velocity field
  allocate(sfd(1:nt)) ! receiver trace
  !allocate(diffop(ny,nx),pmlop(ny,nx))
  allocate(p(ny2,nx2))
  px = 0.
  py = 0.
  ux = 0.
  uy = 0.
  sfd = 0.
  diffop = 0.
  pmlop =0.
  p = 0.

  ! main loop
  do k=2,nt

    ! inject source function
    px(isy,isx) = px(isy,isx) + dt*0.5*fs(k)
    py(isy,isx) = py(isy,isx) + dt*0.5*fs(k)

    do j=1,nx
      do i=1,ny
        ! Update px
        diffop = (ux(i+1,j) - ux(i,j))/dx
        pmlop = qx(i+1,j+1)*px(i+1,j+1)
        px(i+1,j+1) = px(i+1,j+1) - dt*(pmlop + rho*vp(i+1,j+1)*diffop)

        ! Update py
        diffop = (uy(i,j+1) - uy(i,j))/dy
        pmlop = qy(i+1,j+1)*py(i+1,j+1)
        py(i+1,j+1) = py(i+1,j+1) - dt*(pmlop + rho*vp(i+1,j+1)*diffop)

        ! Update ux
        diffop = (px(i+1,j+1) - px(i,j+1) + py(i+1,j+1) - py(i,j+1))/dx
        pmlop = 0.5*(qx(i+1,j+1)+qx(i,j+1))*ux(i,j)
        ux(i,j) = ux(i,j) - dt/rho*(pmlop + diffop)

        ! Update uy
        diffop = (px(i+1,j+1) - px(i+1,j) + py(i+1,j+1) - py(i+1,j))/dy
        pmlop = 0.5*(qy(i+1,j+1)+qy(i+1,j))*uy(i,j)
        uy(i,j) = uy(i,j) - dt/rho*(pmlop + diffop)
      enddo
    enddo

    sfd(k) = px(iry,irx) + py(iry,irx)

    !if (modulo(j,200)==0) then
    !  p = px + py
    !  write (number, *) j !write an internal file to convert an int to a string
    !  open(unit=12, status="replace", form="formatted", file="pressure_"//trim(adjustl(number))//".out")
    !  do k=1,nx2
    !    do l=1,ny2
    !      write(12, *) p(l,k)
    !    enddo
    !  enddo
    !  close(unit=12)
    !endif

  enddo

  ! write output trace in file trace.out
  !open(UNIT=10, STATUS="replace", FORM="formatted", FILE="trace.out")
  !do j=1,nt
  !  write(10, *) sfd(j)
  !enddo
  !close(UNIT=10)

  !open(UNIT=11, STATUS="replace", FORM="formatted", FILE="source.out")
  !do j=1,nt
  !  write(11, *) t(j), fs(j)
  !enddo
  !close(UNIT=11)

  !open(UNIT=12, STATUS="replace", FORM="formatted", FILE="pressure.out")
  !do i=1,nx2
  !  do j=1,ny2
  !  write(12, *) p(j,i)
  !  enddo
  !enddo
  !close(UNIT=12)
end program fdtd2Ds
