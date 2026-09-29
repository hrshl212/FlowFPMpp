program foil_frame_fpm
  ! FPM on a real flow frame: one saved 2-D frame of a flapping-foil simulation
  ! (examples/data/foil_frame_18400.bin, produced by the VICAR3D code and
  ! cropped around the foil and its wake). The foil is represented by the
  ! frame's solid flag as a stair-step body; the influence fields phi_x and
  ! phi_y are solved on the cropped domain with zero outer normal derivative,
  ! and the vortex-induced force densities -(omega x u).grad(phi_i) and their
  ! integrals (per unit span and density) are reported.
  ! Usage: foil_frame_fpm [cpu|gpu [input [output [tolerance]]]]  (one MPI rank)
  ! Input layout (stream, little-endian; read natively, so a little-endian host
  ! such as x86-64 or ARM64 is assumed; the magic string guards the header): character(8) 'FFPM2D01', int32 nx, ny,
  ! real64 xf(0:nx), yf(0:ny), u(nx,ny), v(nx,ny), omega_z(nx,ny), int32 solid(nx,ny).
  ! Output: the same layout as demo_cylinder_vortices (see examples/plot_fields.py).
  use,intrinsic :: iso_fortran_env,only:real64,int32,output_unit
  use mpi
  use fpm_stairstep_mpi_3d
  implicit none
  integer,parameter::nz=2
  character(8)::magic
  character(256)::arg,input,output
  integer(int32)::nx32,ny32
  integer::nx,ny,nproc,ierr,unit,i,j,d,sweeps(2)
  real(real64),allocatable::xf(:),yf(:),xc(:),yc(:),dx(:),dy(:),u2(:,:),v2(:,:),oz2(:,:)
  real(real64),allocatable::phi(:,:,:,:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  real(real64),allocatable::u(:,:,:),v(:,:,:),w(:,:,:),ox(:,:,:),oy(:,:,:),oz(:,:,:),zero(:,:,:)
  real(real64),allocatable::gx(:,:,:),gy(:,:,:),gz(:,:,:),vortex(:,:,:,:),advection(:,:,:)
  integer(int32),allocatable::solid2(:,:)
  integer,allocatable::solid(:,:,:)
  real(real64)::zc(nz),dz(nz),body(3),bounds(6),totals(2,2),correction,imbalance,residual(2),tolerance,t0
  logical::gpu
  gpu=.false.;input='examples/data/foil_frame_18400.bin';output='foil_frame_fpm.bin';tolerance=1e-8_real64
  if(command_argument_count()>0)then;call get_command_argument(1,arg);gpu=trim(arg)=='gpu';endif
  if(command_argument_count()>1)call get_command_argument(2,input)
  if(command_argument_count()>2)call get_command_argument(3,output)
  if(command_argument_count()>3)then;call get_command_argument(4,arg);read(arg,*)tolerance;endif
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=1)error stop 'this example runs on one MPI rank'
  open(newunit=unit,file=trim(input),access='stream',form='unformatted',status='old')
  read(unit)magic,nx32,ny32
  if(magic/='FFPM2D01')error stop 'unexpected input format'
  if(nx32<2.or.ny32<2.or.nx32>100000.or.ny32>100000)error stop 'implausible grid size (byte order?)'
  nx=nx32;ny=ny32
  allocate(xf(0:nx),yf(0:ny),u2(nx,ny),v2(nx,ny),oz2(nx,ny),solid2(nx,ny))
  read(unit)xf,yf,u2,v2,oz2,solid2
  close(unit)
  ! Cell centres with one mirrored halo centre per side, and widths.
  allocate(xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny))
  dx=xf(1:nx)-xf(0:nx-1);dy=yf(1:ny)-yf(0:ny-1)
  xc(1:nx)=0.5_real64*(xf(0:nx-1)+xf(1:nx));yc(1:ny)=0.5_real64*(yf(0:ny-1)+yf(1:ny))
  xc(0)=xc(1)-dx(1);xc(nx+1)=xc(nx)+dx(nx);yc(0)=yc(1)-dy(1);yc(ny+1)=yc(ny)+dy(ny)
  ! Two z cells of unit total span with zero z-face derivatives make the case 2-D.
  dz=0.5_real64;zc=[0.25_real64,0.75_real64]
  allocate(phi(0:nx+1,0:ny+1,nz,2),rhs(nx,ny,nz),solid(0:nx+1,0:ny+1,nz))
  allocate(xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny))
  allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz),zero(nx,ny,nz))
  allocate(gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz),vortex(nx,ny,nz,2),advection(nx,ny,nz))
  solid=0;w=0;ox=0;oy=0;zero=0
  do d=1,nz
    solid(1:nx,1:ny,d)=solid2;u(:,:,d)=u2;v(:,:,d)=v2;oz(:,:,d)=oz2
  enddo
  where(solid(1:nx,1:ny,:)/=0)
    u=0;v=0;oz=0
  endwhere
  bounds=[-huge(1.0_real64),huge(1.0_real64),-huge(1.0_real64),huge(1.0_real64),-huge(1.0_real64),huge(1.0_real64)]
  write(output_unit,'(a,i0,a,i0,a,i0,a)')'frame: ',nx,' x ',ny,' cells, ',count(solid2/=0),' solid'
  write(output_unit,'(a,4(f10.4,1x))')'extent x0 x1 y0 y1: ',xf(0),xf(nx),yf(0),yf(ny)
  rhs=0;t0=MPI_Wtime()
  do d=1,2
    body=0;body(d)=1;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0;phi(:,:,:,d)=0
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD,correction,imbalance)
    call solve_stairstep_poisson_3d(phi(:,:,:,d),rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD, &
         tolerance,2000000,50,sweeps(d),residual(d),use_openacc=gpu,method='pcg')
    call compute_stairstep_fpm_terms_3d(phi(:,:,:,d),u,v,w,ox,oy,oz,zero,zero,zero,solid,body,xc,yc,zc,dx,dy,dz, &
         xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vortex(:,:,:,d),advection,totals(:,d), &
         nx,ny,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD,use_openacc=gpu)
    write(output_unit,'(a,i0,a,i0,a,es10.3,a,es10.3)')'direction ',d,': PCG iterations ',sweeps(d), &
         ', residual ',residual(d),', compatibility correction ',correction
  enddo
  write(output_unit,'(a,f10.3)')'solve and terms seconds: ',MPI_Wtime()-t0
  ! The z span is 1, so the integrals are forces per unit span (and density).
  write(output_unit,'(a,2(es16.8,1x))')'vortex-induced force (x, y): ',totals(1,1),totals(1,2)
  open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  write(unit)int(nx,int32),int(ny,int32),xc(1:nx),yc(1:ny),oz(:,:,1),vortex(:,:,1,1),vortex(:,:,1,2), &
       phi(1:nx,1:ny,1,1),phi(1:nx,1:ny,1,2),real(solid(1:nx,1:ny,1),real64)
  close(unit)
  write(output_unit,'(a)')'wrote '//trim(output)
  call MPI_Finalize(ierr)
end program foil_frame_fpm
