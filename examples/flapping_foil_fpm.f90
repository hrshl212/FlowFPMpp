program flapping_foil_fpm
  ! Force partitioning of one frame of a 2-D passive-flapping foil simulation
  ! (examples/data/flapping_foil_216000.bin: step 216000, t = 172.8, of a VICAR3D
  ! run of a heaving and passively pitching foil with an angle limiter, Re = 1e4).
  ! The foil is the frame's solid flag, treated as a stair-step body. The
  ! influence fields phi_x and phi_y are solved on the full domain with zero outer
  ! normal derivative (PCG), and the four terms of the force decomposition
  ! (F_KINE, F_VIF, F_VIS, F_BOUND; see src/fpm_partition_2d.f90) are reported
  ! per unit span and unit density, in the simulation's units (U = 1, chord = 1).
  ! Usage: flapping_foil_fpm [cpu|gpu [input [output [tolerance]]]]  (one MPI rank)
  ! Input (Fortran stream, little-endian, read natively; the magic string and a
  ! size check guard the header):
  !   character(8) 'FFPM2D02'; int32 nx, ny; real64 xf(0:nx), yf(0:ny)
  !   real64 u(nx,ny), v(nx,ny); int8 solid(nx,ny); real64 dt
  !   real64 ring(2*nx+2*ny, 2, 2)   outer-boundary cells (west column, east
  !        column, south row, north row) x (u, v) x (t - dt, t + dt)
  !   real64 kin(6)   body velocity (k1 - k3 y, k2 + k3 x) and acceleration
  !        (k4 - k6 y - k3^2 x, k5 + k6 x - k3^2 y); real64 time
  ! Output: the layout of demo_cylinder_vortices, with the F_VIF densities in
  ! place of the force densities (see examples/plot_fields.py).
  use,intrinsic :: iso_fortran_env,only:real64,int8,int32,output_unit
  use mpi
  use fpm_stairstep_mpi_3d
  use fpm_partition_2d
  implicit none
  integer,parameter::nz=2
  real(real64),parameter::re=1.0e4_real64
  character(8)::magic
  character(256)::arg,input,output
  integer(int32)::nx32,ny32
  integer::nx,ny,nproc,ierr,unit,d,k,sweeps(2)
  real(real64),allocatable::xf(:),yf(:),xc(:),yc(:),dx(:),dy(:),u2(:,:),v2(:,:),ring(:,:,:)
  real(real64),allocatable::phi(:,:,:,:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  real(real64),allocatable::vif(:,:,:),omega(:,:)
  integer(int8),allocatable::solid1(:,:)
  integer,allocatable::solid(:,:,:)
  logical,allocatable::lsolid(:,:)
  real(real64)::zc(nz),dz(nz),body(3),kin(6),dt,time,terms(4,2),correction,imbalance,residual(2),tolerance,t0
  logical::gpu
  gpu=.false.;input='examples/data/flapping_foil_216000.bin';output='flapping_foil_fpm.bin';tolerance=1e-10_real64
  if(command_argument_count()>0)then;call get_command_argument(1,arg);gpu=trim(arg)=='gpu';endif
  if(command_argument_count()>1)call get_command_argument(2,input)
  if(command_argument_count()>2)call get_command_argument(3,output)
  if(command_argument_count()>3)then;call get_command_argument(4,arg);read(arg,*)tolerance;endif
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=1)error stop 'this example runs on one MPI rank'
  open(newunit=unit,file=trim(input),access='stream',form='unformatted',status='old')
  read(unit)magic,nx32,ny32
  if(magic/='FFPM2D02')error stop 'unexpected input format'
  if(nx32<2.or.ny32<2.or.nx32>100000.or.ny32>100000)error stop 'implausible grid size (byte order?)'
  nx=nx32;ny=ny32
  allocate(xf(0:nx),yf(0:ny),u2(nx,ny),v2(nx,ny),solid1(nx,ny),ring(2*nx+2*ny,2,2))
  read(unit)xf,yf,u2,v2,solid1,dt,ring,kin,time
  close(unit)
  allocate(xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny),lsolid(nx,ny))
  lsolid=solid1/=0
  dx=xf(1:nx)-xf(0:nx-1);dy=yf(1:ny)-yf(0:ny-1)
  xc(1:nx)=0.5_real64*(xf(0:nx-1)+xf(1:nx));yc(1:ny)=0.5_real64*(yf(0:ny-1)+yf(1:ny))
  xc(0)=xc(1)-dx(1);xc(nx+1)=xc(nx)+dx(nx);yc(0)=yc(1)-dy(1);yc(ny+1)=yc(ny)+dy(ny)
  ! Two z cells of unit total span with zero z-face derivatives make the case 2-D.
  dz=0.5_real64;zc=[0.25_real64,0.75_real64]
  allocate(phi(0:nx+1,0:ny+1,nz,2),rhs(nx,ny,nz),solid(0:nx+1,0:ny+1,nz))
  allocate(xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny),vif(nx,ny,2),omega(nx,ny))
  solid=0
  do k=1,nz;solid(1:nx,1:ny,k)=merge(1,0,lsolid);enddo
  write(output_unit,'(a,i0,a,i0,a,i0,a,f8.3)')'frame: ',nx,' x ',ny,' cells, ',count(lsolid),' solid, t = ',time
  write(output_unit,'(a,3es12.4)')'body velocity U_B = (k1 - k3 y, k2 + k3 x), k1..k3: ',kin(1:3)
  write(output_unit,'(a,3es12.4)')'body acceleration terms k4..k6: ',kin(4:6)
  rhs=0;t0=MPI_Wtime()
  do d=1,2
    body=0;body(d)=1;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0;phi(:,:,:,d)=0
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD,correction,imbalance)
    call solve_stairstep_poisson_3d(phi(:,:,:,d),rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD, &
         tolerance,2000000,50,sweeps(d),residual(d),use_openacc=gpu,method='pcg')
    call fpm_partition_terms_2d(phi(1:nx,1:ny,1,d),body(1:2),u2,v2,lsolid,xf,yf,kin,ring,dt,1/re, &
         nx,ny,terms(:,d),vif(:,:,d),omega)
    write(output_unit,'(a,i0,a,i0,a,es10.3)')'direction ',d,': PCG iterations ',sweeps(d),', residual ',residual(d)
  enddo
  write(output_unit,'(a,f10.3)')'solve and partition seconds: ',MPI_Wtime()-t0
  write(output_unit,'(a)')'force partition per unit span and density (x = drag direction, y = lift):'
  write(output_unit,'(a,2es14.5)')'  F_KINE  (body acceleration)      ',terms(1,:)
  write(output_unit,'(a,2es14.5)')'  F_VIF   (vortex-induced)         ',terms(2,:)
  write(output_unit,'(a,2es14.5)')'  F_VIS   (viscous)                ',terms(3,:)
  write(output_unit,'(a,2es14.5)')'  F_BOUND (outer boundary)         ',terms(4,:)
  write(output_unit,'(a,2es14.5)')'  sum = pressure force             ',sum(terms,dim=1)
  open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  write(unit)int(nx,int32),int(ny,int32),xc(1:nx),yc(1:ny),omega,vif(:,:,1),vif(:,:,2), &
       phi(1:nx,1:ny,1,1),phi(1:nx,1:ny,1,2),merge(1.0_real64,0.0_real64,lsolid)
  close(unit)
  write(output_unit,'(a)')'wrote '//trim(output)
  call MPI_Finalize(ierr)
end program flapping_foil_fpm
