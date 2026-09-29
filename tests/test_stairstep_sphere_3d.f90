program test_stairstep_sphere_3d
  ! Grid-convergence study of the stair-step influence field of a sphere.
  ! For a sphere of radius a the x influence potential with n.grad(phi)=n_x on
  ! the body and decay at infinity is the dipole phi=-a^3 x/(2 r^3); its exact
  ! gradient is imposed on the outer faces of the box [-1.5,1.5]^3. The test
  ! reports the fluid-volume RMS and far-field (r>2a) maximum error of phi,
  ! and the added-mass coefficient C_m=-(1/V) sum(phi n_x dA) over body faces
  ! (exact value 1/2), on uniform grids N=12,24,48 (and 96 with a second
  ! argument "extended"). Stair-stepping makes the method first order at best;
  ! the test requires monotone decrease of all three errors and an observed
  ! order above 0.7 for the field errors on the finest refinement.
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use mpi
  use fpm_test_grid
  use fpm_stairstep_mpi_3d
  implicit none
  real(real64),parameter::a=0.5_real64,half=1.5_real64
  integer,parameter::levels(4)=[12,24,48,96]
  real(real64)::rms(4),far(4),cm(4)
  integer::rank,ierr,m,sweeps(4),nl
  logical::use_gpu,passed
  character(16)::arg
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  use_gpu=.false.
  if(command_argument_count()>0)then;call get_command_argument(1,arg);use_gpu=trim(arg)=='openacc';endif
  nl=3
  if(command_argument_count()>1)then;call get_command_argument(2,arg);if(trim(arg)=='extended')nl=4;endif
  do m=1,nl
    call run_level(levels(m),rms(m),far(m),cm(m),sweeps(m))
  enddo
  if(rank==0)then
    write(output_unit,'(a)')'   N  sweeps   RMS error   far-field Linf   C_m (exact 0.5)'
    do m=1,nl
      write(output_unit,'(i4,i8,2(es13.4),f12.5)')levels(m),sweeps(m),rms(m),far(m),cm(m)
    enddo
    write(output_unit,'(a,3f8.3)')'observed RMS orders: ',log(rms(1:nl-1)/rms(2:nl))/log(2.0_real64)
    write(output_unit,'(a,3f8.3)')'observed far-field orders: ',log(far(1:nl-1)/far(2:nl))/log(2.0_real64)
    write(output_unit,'(a,3f8.3)')'observed C_m orders: ', &
         log(abs(cm(1:nl-1)-0.5_real64)/abs(cm(2:nl)-0.5_real64))/log(2.0_real64)
    flush(output_unit)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  passed=all(rms(2:nl)<rms(1:nl-1)).and.all(far(2:nl)<far(1:nl-1))
  passed=passed.and.all(ieee_is_finite([rms(1:nl),far(1:nl),cm(1:nl)]))
  passed=passed.and.all(abs(cm(2:nl)-0.5_real64)<abs(cm(1:nl-1)-0.5_real64))
  passed=passed.and.log(rms(nl-1)/rms(nl))/log(2.0_real64)>0.7_real64
  passed=passed.and.log(far(nl-1)/far(nl))/log(2.0_real64)>0.7_real64
  if(.not.passed)then
    if(rank==0)write(error_unit,'(a)')'FAIL: stair-step sphere convergence'
    flush(output_unit);flush(error_unit)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif
  if(rank==0)write(output_unit,'(a,l1,a)')'PASS: stair-step sphere convergence (openacc=',use_gpu,')'
  call MPI_Finalize(ierr)
contains
  subroutine run_level(n,rms,far,cm,sweeps)
    integer,intent(in)::n
    real(real64),intent(out)::rms,far,cm
    integer,intent(out)::sweeps
    real(real64),allocatable::gf(:),gc(:),gw(:),zf(:),zc(:),dz(:),xc(:),yc(:),dx(:),dy(:)
    real(real64),allocatable::phi(:,:,:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
    integer,allocatable::solid(:,:,:)
    real(real64)::body(3),correction,imbalance,residual,local(4),global(4),mean,e,vol,r
    integer::cart,nx,ny,ioff,joff,west,east,south,north,i,j,k
    allocate(gf(0:n),gc(n),gw(n),zf(0:n),zc(n),dz(n))
    call stretched_axis(n,-half,half,0.0_real64,gf,gc,gw)
    zf=gf;zc=gc;dz=gw
    call cartesian_layout(MPI_COMM_WORLD,n,n,cart,nx,ny,ioff,joff,west,east,south,north)
    allocate(xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny),phi(0:nx+1,0:ny+1,n),rhs(nx,ny,n),solid(0:nx+1,0:ny+1,n))
    allocate(xlo(ny,n),xhi(ny,n),ylo(nx,n),yhi(nx,n),zlo(nx,ny),zhi(nx,ny))
    call local_axis(n,gc,gw,ioff,nx,xc,dx)
    call local_axis(n,gc,gw,joff,ny,yc,dy)
    solid=0
    do k=1,n;do j=1,ny;do i=1,nx
      if(xc(i)**2+yc(j)**2+zc(k)**2<a*a)solid(i,j,k)=1
    enddo;enddo;enddo
    call exchange_solid_halo_3d(solid,nx,ny,n,west,east,south,north,cart)
    do k=1,n;do j=1,ny
      xlo(j,k)=dipole_grad(1,-half,yc(j),zc(k));xhi(j,k)=dipole_grad(1,half,yc(j),zc(k))
    enddo;enddo
    do k=1,n;do i=1,nx
      ylo(i,k)=dipole_grad(2,xc(i),-half,zc(k));yhi(i,k)=dipole_grad(2,xc(i),half,zc(k))
    enddo;enddo
    do j=1,ny;do i=1,nx
      zlo(i,j)=dipole_grad(3,xc(i),yc(j),-half);zhi(i,j)=dipole_grad(3,xc(i),yc(j),half)
    enddo;enddo
    rhs=0;phi=0;body=[1.0_real64,0.0_real64,0.0_real64]
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,n,west,east,south,north,cart,correction,imbalance)
    call solve_stairstep_poisson_3d(phi,rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,n,west,east,south,north,cart,1.0e-10_real64,2000000,100,sweeps,residual,use_openacc=use_gpu)
    ! Error relative to the dipole after removing the fluid mean of the difference.
    local=0
    do k=1,n;do j=1,ny;do i=1,nx
      if(solid(i,j,k)/=0)cycle
      vol=dx(i)*dy(j)*dz(k)
      local(1:2)=local(1:2)+[vol*(phi(i,j,k)-dipole(xc(i),yc(j),zc(k))),vol]
    enddo;enddo;enddo
    call MPI_Allreduce(local,global,2,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
    mean=global(1)/global(2)
    local=0
    do k=1,n;do j=1,ny;do i=1,nx
      if(solid(i,j,k)/=0)cycle
      vol=dx(i)*dy(j)*dz(k)
      e=phi(i,j,k)-dipole(xc(i),yc(j),zc(k))-mean
      r=sqrt(xc(i)**2+yc(j)**2+zc(k)**2)
      local(1)=local(1)+vol*e*e
      if(r>2*a)local(2)=max(local(2),abs(e))
      ! Body x-faces: the outward body normal is -x on the fluid cell's east
      ! side. phi is extrapolated to the face with the imposed gradient (1).
      if(i<nx.or.east/=MPI_PROC_NULL)then
        if(solid(i+1,j,k)/=0)local(3)=local(3)-(phi(i,j,k)+0.5_real64*dx(i))*dy(j)*dz(k)
      endif
      if(i>1.or.west/=MPI_PROC_NULL)then
        if(solid(i-1,j,k)/=0)local(3)=local(3)+(phi(i,j,k)-0.5_real64*dx(i))*dy(j)*dz(k)
      endif
    enddo;enddo;enddo
    call MPI_Allreduce(local(1),global(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
    call MPI_Allreduce(local(2),global(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,cart,ierr)
    call MPI_Allreduce(local(3),global(3),1,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
    rms=sqrt(global(1)/volume_fluid(solid,dx,dy,dz,nx,ny,n,cart));far=global(2)
    cm=-global(3)/(4*acos(-1.0_real64)*a**3/3)
    call MPI_Comm_free(cart,ierr)
  end subroutine run_level

  real(real64) function volume_fluid(solid,dx,dy,dz,nx,ny,nz,comm)
    integer,intent(in)::nx,ny,nz,comm,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::dx(nx),dy(ny),dz(nz)
    real(real64)::local
    integer::i,j,k,ierr
    local=0
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)==0)local=local+dx(i)*dy(j)*dz(k)
    enddo;enddo;enddo
    call MPI_Allreduce(local,volume_fluid,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  end function volume_fluid

  pure real(real64) function dipole(x,y,z)
    real(real64),intent(in)::x,y,z
    dipole=-a**3*x/(2*sqrt(x*x+y*y+z*z)**3)
  end function dipole

  pure real(real64) function dipole_grad(axis,x,y,z)
    integer,intent(in)::axis
    real(real64),intent(in)::x,y,z
    real(real64)::r2,c(3)
    r2=x*x+y*y+z*z;c=[x,y,z]
    dipole_grad=-0.5_real64*a**3*(merge(1.0_real64,0.0_real64,axis==1)/sqrt(r2)**3-3*x*c(axis)/sqrt(r2)**5)
  end function dipole_grad
end program test_stairstep_sphere_3d
