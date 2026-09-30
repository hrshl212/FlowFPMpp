program test_fpm_partition_2d
  ! Checks of the force-partitioning terms (src/fpm_partition_2d.f90):
  ! 1. solid-body rotation u = W(-y, x): div[(u.grad)u] = -2 W^2 exactly in
  !    interior cells (linear face fluxes are exact);
  ! 2. uniform flow with no body: every term vanishes exactly;
  ! 3. a stair-step cylinder of radius R accelerating at a = (1, 0) through
  !    fluid at rest: F_KINE = -m_a with the 2-D added mass pi R^2 (checked to
  !    within 20%, allowing for the stair-step body and the finite box), and the other
  !    terms vanish. The stair-step body overestimates the added mass (it
  !    converges from above, as in the sphere study), so 20% is allowed at R = 10h.
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use mpi
  use fpm_test_grid
  use fpm_stairstep_mpi_3d
  use fpm_partition_2d
  implicit none
  integer,parameter::n=160,nz=2
  real(real64),parameter::half=2.0_real64,radius=0.25_real64,w0=0.7_real64
  real(real64)::xf(0:n),xcg(n),dxg(n),xc(0:n+1),dx(n),zc(nz),dz(nz)
  real(real64),allocatable::u(:,:),v(:,:),phi(:,:),ring(:,:,:),dens(:,:),omega(:,:)
  real(real64),allocatable::phi3(:,:,:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  integer,allocatable::solid3(:,:,:)
  logical,allocatable::solid(:,:)
  real(real64)::terms(4),kin(6),err1,err2,cm,correction,imbalance,residual,body(3)
  integer::i,j,ierr,sweeps
  logical::passed
  call MPI_Init(ierr)
  call stretched_axis(n,-half,half,0.0_real64,xf,xcg,dxg)
  call local_axis(n,xcg,dxg,0,n,xc,dx)
  allocate(u(n,n),v(n,n),phi(n,n),ring(4*n,2,2),dens(n,n),omega(n,n),solid(n,n))
  ! 1. Solid-body rotation, no body, phi = 1.
  solid=.false.;phi=1;kin=0
  do j=1,n;do i=1,n;u(i,j)=-w0*xc(j);v(i,j)=w0*xc(i);enddo;enddo
  call fill_ring()
  call fpm_partition_terms_2d(phi,[1.0_real64,0.0_real64],u,v,solid,xf,xf,kin,ring,1.0_real64,1.0e-3_real64, &
       n,n,terms,dens,omega)
  err1=maxval(abs(dens(3:n-2,3:n-2)+2*w0**2))
  err1=max(err1,maxval(abs(omega(3:n-2,3:n-2)-2*w0)))
  ! 2. Uniform flow, no body.
  u=1;v=0;call fill_ring()
  do j=1,n;phi(:,j)=xc(1:n)-xc(j);enddo
  call fpm_partition_terms_2d(phi,[1.0_real64,0.0_real64],u,v,solid,xf,xf,kin,ring,1.0_real64,1.0e-3_real64, &
       n,n,terms,dens,omega)
  err2=maxval(abs(terms))
  ! 3. Accelerating cylinder in fluid at rest.
  allocate(phi3(0:n+1,0:n+1,nz),rhs(n,n,nz),solid3(0:n+1,0:n+1,nz))
  allocate(xlo(n,nz),xhi(n,nz),ylo(n,nz),yhi(n,nz),zlo(n,n),zhi(n,n))
  dz=0.5_real64;zc=[0.25_real64,0.75_real64]
  do j=1,n;do i=1,n;solid(i,j)=xc(i)**2+xc(j)**2<radius**2;enddo;enddo
  solid3=0;do i=1,nz;solid3(1:n,1:n,i)=merge(1,0,solid);enddo
  rhs=0;phi3=0;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0;body=[1.0_real64,0.0_real64,0.0_real64]
  call enforce_stairstep_compatibility_3d(rhs,solid3,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dx,dz, &
       n,n,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD,correction,imbalance)
  call solve_stairstep_poisson_3d(phi3,rhs,solid3,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,xc,zc,dx,dx,dz, &
       n,n,nz,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_COMM_WORLD, &
       1.0e-10_real64,200000,50,sweeps,residual,use_openacc=.false.)
  u=0;v=0;call fill_ring();kin=[0.0_real64,0.0_real64,0.0_real64,1.0_real64,0.0_real64,0.0_real64]
  call fpm_partition_terms_2d(phi3(1:n,1:n,1),[1.0_real64,0.0_real64],u,v,solid,xf,xf,kin,ring,1.0_real64, &
       1.0e-3_real64,n,n,terms,dens,omega)
  cm=-terms(1)/(acos(-1.0_real64)*radius**2)
  write(output_unit,'(a,es12.4)')'rotation: max |div N + 2W^2|, |omega - 2W| (interior) = ',err1
  write(output_unit,'(a,es12.4)')'uniform flow: max |term| = ',err2
  write(output_unit,'(a,i0,a,es10.3)')'cylinder: PCG iterations ',sweeps,', residual ',residual
  write(output_unit,'(a,4es12.4)')'cylinder terms KINE VIF VIS BOUND = ',terms
  write(output_unit,'(a,f8.4)')'cylinder added-mass coefficient -F_KINE/(pi R^2) = ',cm
  passed=err1<1e-10.and.err2==0.and.abs(cm-1)<0.2_real64.and.maxval(abs(terms(2:4)))==0
  passed=passed.and.all(ieee_is_finite([err1,err2,cm,terms]))
  if(.not.passed)then
    write(error_unit,'(a)')'FAIL: force-partitioning terms';call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif
  write(output_unit,'(a)')'PASS: force-partitioning terms'
  call MPI_Finalize(ierr)
contains
  subroutine fill_ring()
    integer::t
    do t=1,2
      ring(:,1,t)=[u(1,:),u(n,:),u(:,1),u(:,n)]
      ring(:,2,t)=[v(1,:),v(n,:),v(:,1),v(:,n)]
    enddo
  end subroutine fill_ring
end program test_fpm_partition_2d
