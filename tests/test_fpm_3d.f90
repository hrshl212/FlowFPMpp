program test_fpm_3d
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use fpm_laplace_3d
  use fpm_terms_3d
  implicit none
  integer,parameter::nx=12,ny=10,nz=8
  real(real64),parameter::tol=1.0e-10_real64
  real(real64),allocatable::phi(:,:,:),rhs(:,:,:),exact(:,:,:),u(:,:,:),v(:,:,:),w(:,:,:)
  real(real64),allocatable::ox(:,:,:),oy(:,:,:),oz(:,:,:),ax(:,:,:),ay(:,:,:),az(:,:,:)
  real(real64),allocatable::gx(:,:,:),gy(:,:,:),gz(:,:,:),vf(:,:,:),af(:,:,:)
  real(real64),allocatable::xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  integer,allocatable::fluid(:,:,:)
  real(real64)::dx,dy,dz,x,y,z,mean_exact,correction,imbalance,residual,totals(2)
  real(real64)::solution_error,gradient_error,field_error,expected_v,expected_a
  integer::i,j,k,iterations
  logical::passed
  dx=1.0_real64/nx;dy=1.0_real64/ny;dz=1.0_real64/nz
  allocate(phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz),exact(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz), &
       ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz),ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz), &
       gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz),vf(nx,ny,nz),af(nx,ny,nz),fluid(nx,ny,nz), &
       xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny))
  phi=0;rhs=1;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0
  call enforce_neumann_compatibility_3d(rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz,correction,imbalance)
  do k=1,nz;z=(k-0.5_real64)*dz;do j=1,ny;y=(j-0.5_real64)*dy;do i=1,nx;x=(i-0.5_real64)*dx
    exact(i,j,k)=((x-.5_real64)**2+(y-.5_real64)**2+(z-.5_real64)**2)/6.0_real64
  enddo;enddo;enddo
  mean_exact=sum(exact)/real(nx*ny*nz,real64);exact=exact-mean_exact
  call solve_neumann_poisson_3d(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz,tol,50000,25,iterations,residual)
  u=1;v=2;w=3;ox=0.5_real64;oy=-0.25_real64;oz=0.75_real64
  ax=0.2_real64;ay=-0.1_real64;az=0.3_real64;fluid=1
  call compute_fpm_terms_uniform_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,gx,gy,gz,vf,af,totals,nx,ny,nz,dx,dy,dz)
  solution_error=maxval(abs(phi(1:nx,1:ny,1:nz)-exact));gradient_error=0;field_error=0
  do k=1,nz;z=(k-.5_real64)*dz;do j=1,ny;y=(j-.5_real64)*dy;do i=1,nx;x=(i-.5_real64)*dx
    gradient_error=max(gradient_error,abs(gx(i,j,k)-(x-.5_real64)/3),abs(gy(i,j,k)-(y-.5_real64)/3), &
         abs(gz(i,j,k)-(z-.5_real64)/3))
    expected_v=-((-2.25_real64)*((x-.5_real64)/3)-.75_real64*((y-.5_real64)/3)+1.25_real64*((z-.5_real64)/3))
    expected_a=-(.2_real64*((x-.5_real64)/3)-.1_real64*((y-.5_real64)/3)+.3_real64*((z-.5_real64)/3))
    field_error=max(field_error,abs(vf(i,j,k)-expected_v),abs(af(i,j,k)-expected_a))
  enddo;enddo;enddo
  write(output_unit,'(a,es14.6)')'compatibility correction: ',correction
  write(output_unit,'(a,es14.6)')'compatibility imbalance: ',imbalance
  write(output_unit,'(a,i0,a,es14.6)')'iterations: ',iterations,', residual: ',residual
  write(output_unit,'(a,es14.6)')'solution error (Linf): ',solution_error
  write(output_unit,'(a,es14.6)')'gradient error (Linf): ',gradient_error
  write(output_unit,'(a,es14.6)')'FPM field error (Linf): ',field_error
  write(output_unit,'(a,2(es14.6,1x))')'integrated terms: ',totals
  passed=iterations<50000.and.residual<=tol.and.abs(correction-1.0_real64/6)<1e-12.and.abs(imbalance)<1e-12
  passed=passed.and.solution_error<1e-9.and.gradient_error<1e-9.and.field_error<1e-9.and.maxval(abs(totals))<1e-12
  if(.not.passed)then;write(error_unit,'(a)')'FAIL: manufactured 3-D FPM test';error stop 1;endif
  write(output_unit,'(a)')'PASS: manufactured 3-D Laplace/FPM test'
end program test_fpm_3d
