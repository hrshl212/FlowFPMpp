program demo_cylinder_vortices
  ! Demonstration: FPM vortex-induced force (divergence form, see
  ! src/fpm_partition_2d.f90) around a stair-step circular
  ! cylinder (radius 0.5, centre at the origin) in a synthetic flow made of a
  ! uniform stream plus a staggered row of Lamb-Oseen vortices. The flow is a
  ! prescribed kinematic input, not a Navier-Stokes solution. The case is 2-D:
  ! two z cells with zero z-face derivatives. The cylinder is stationary and the
  ! flow steady and treated as inviscid, so F_KINE and F_VIS vanish. The
  ! prescribed flow does not satisfy no-slip, so the body faces are treated as
  ! slip faces (slip_body) to avoid a spurious wall layer. Since the
  ! flow is not a Navier-Stokes solution, the terms are not balanced by a
  ! pressure force.
  ! Usage: demo_cylinder_vortices [cpu|gpu [gnx gny [output]]]  (one MPI rank)
  ! Output (stream, native endian): nx, ny (int32), xc(nx), yc(ny), then
  ! omega_z, F_VIF density x, F_VIF density y, phi_x, phi_y, solid, each nx*ny doubles.
  use,intrinsic :: iso_fortran_env,only:real64,int32,output_unit
  use mpi
  use fpm_test_grid
  use fpm_stairstep_mpi_3d
  use fpm_partition_2d
  implicit none
  integer,parameter::nz=2,nvort=6
  real(real64),parameter::radius=0.5_real64,core=0.22_real64
  real(real64)::vx(nvort),vy(nvort),circ(nvort)
  real(real64),allocatable::xf(:),xcg(:),dxg(:),yf(:),ycg(:),dyg(:),xc(:),yc(:),dx(:),dy(:)
  real(real64),allocatable::phi(:,:,:,:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  real(real64),allocatable::u(:,:,:),v(:,:,:),w(:,:,:),ox(:,:,:),oy(:,:,:),oz(:,:,:),zero(:,:,:)
  real(real64),allocatable::vortex(:,:,:,:)
  integer,allocatable::solid(:,:,:)
  real(real64)::zf(0:nz),zc(nz),dz(nz),body(3),totals(4,2),correction,imbalance,residual,t0,r2,q
  real(real64),allocatable::ring(:,:,:),omega2(:,:)
  real(real64),parameter::kin0(6)=0.0_real64
  logical,allocatable::lsolid(:,:)
  integer::gnx,gny,nproc,ierr,cart,nx,ny,ioff,joff,west,east,south,north,i,j,m,d,sweeps(2),unit
  logical::gpu
  character(256)::arg,output
  gpu=.false.;gnx=360;gny=240;output='demo_cylinder_vortices.bin'
  if(command_argument_count()>0)then;call get_command_argument(1,arg);gpu=trim(arg)=='gpu';endif
  if(command_argument_count()>=3)then
    call get_command_argument(2,arg);read(arg,*)gnx
    call get_command_argument(3,arg);read(arg,*)gny
  endif
  if(command_argument_count()>=4)call get_command_argument(4,output)
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=1)error stop 'the demo writes one global file; run it on one MPI rank'
  call cartesian_layout(MPI_COMM_WORLD,gnx,gny,cart,nx,ny,ioff,joff,west,east,south,north)
  allocate(xf(0:gnx),xcg(gnx),dxg(gnx),yf(0:gny),ycg(gny),dyg(gny),xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny))
  call stretched_axis(gnx,-3.0_real64,6.0_real64,0.0_real64,xf,xcg,dxg)
  call stretched_axis(gny,-3.0_real64,3.0_real64,0.0_real64,yf,ycg,dyg)
  call stretched_axis(nz,0.0_real64,2*dxg(1),0.0_real64,zf,zc,dz)
  call local_axis(gnx,xcg,dxg,ioff,nx,xc,dx)
  call local_axis(gny,ycg,dyg,joff,ny,yc,dy)
  allocate(phi(0:nx+1,0:ny+1,nz,2),rhs(nx,ny,nz),solid(0:nx+1,0:ny+1,nz))
  allocate(xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny))
  allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz),zero(nx,ny,nz))
  allocate(vortex(nx,ny,nz,2),ring(2*nx+2*ny,2,2),omega2(nx,ny),lsolid(nx,ny))
  ! Staggered vortex row in the near wake, alternating sign.
  do m=1,nvort
    vx(m)=0.9_real64+0.75_real64*(m-1);vy(m)=merge(0.3_real64,-0.3_real64,mod(m,2)==1)
    circ(m)=merge(-1.6_real64,1.6_real64,mod(m,2)==1)*(1-0.08_real64*(m-1))
  enddo
  solid=0;u=1;v=0;w=0;ox=0;oy=0;oz=0;zero=0
  do j=1,ny;do i=1,nx
    if(xc(i)**2+yc(j)**2<radius**2)solid(i,j,:)=1
    do m=1,nvort
      r2=(xc(i)-vx(m))**2+(yc(j)-vy(m))**2
      q=circ(m)/(2*acos(-1.0_real64))*(1-exp(-r2/core**2))/max(r2,1e-30_real64)
      u(i,j,:)=u(i,j,:)-q*(yc(j)-vy(m));v(i,j,:)=v(i,j,:)+q*(xc(i)-vx(m))
      oz(i,j,:)=oz(i,j,:)+circ(m)/(acos(-1.0_real64)*core**2)*exp(-r2/core**2)
    enddo
  enddo;enddo
  call exchange_solid_halo_3d(solid,nx,ny,nz,west,east,south,north,cart)
  lsolid=solid(1:nx,1:ny,1)/=0
  ! Steady flow: identical outer-boundary snapshots, so du/dt = 0.
  do d=1,2
    ring(:,1,d)=[u(1,:,1),u(nx,:,1),u(:,1,1),u(:,ny,1)]
    ring(:,2,d)=[v(1,:,1),v(nx,:,1),v(:,1,1),v(:,ny,1)]
  enddo
  rhs=0;t0=MPI_Wtime()
  do d=1,2
    body=0;body(d)=1;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0;phi(:,:,:,d)=0
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,correction,imbalance)
    call solve_stairstep_poisson_3d(phi(:,:,:,d),rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,1.0e-9_real64,5000000,500,sweeps(d),residual,use_openacc=gpu)
    call fpm_partition_terms_2d(phi(1:nx,1:ny,1,d),body(1:2),u(:,:,1),v(:,:,1),lsolid,xf,yf,kin0, &
         ring,1.0_real64,0.0_real64,nx,ny,totals(:,d),vortex(:,:,1,d),omega2,slip_body=.true.)
    write(output_unit,'(a,i0,a,i0,a,es11.3,a,es11.3)')'direction ',d,': sweeps ',sweeps(d), &
         ', residual ',residual,', compatibility correction ',correction
  enddo
  write(output_unit,'(a,2(es14.6,1x))')'vortex-induced force F_VIF per unit span (x, y): ',totals(2,1),totals(2,2)
  write(output_unit,'(a,f10.3)')'solve and terms seconds: ',MPI_Wtime()-t0
  open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  write(unit)int(nx,int32),int(ny,int32),xc(1:nx),yc(1:ny),oz(:,:,1),vortex(:,:,1,1),vortex(:,:,1,2), &
       phi(1:nx,1:ny,1,1),phi(1:nx,1:ny,1,2),real(solid(1:nx,1:ny,1),real64)
  close(unit)
  write(output_unit,'(a)')'wrote '//trim(output)
  call MPI_Comm_free(cart,ierr)
  call MPI_Finalize(ierr)
end program demo_cylinder_vortices
