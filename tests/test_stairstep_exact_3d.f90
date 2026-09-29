program test_stairstep_exact_3d
  ! Stair-step solver checks on a stretched grid with several blocked bodies
  ! (one touching the outer boundary), for any rank count whose 2-D layout
  ! divides 24x20. Pass "openacc" to run the distributed solves on the device.
  ! 1. Exactness: with body and outer face gradients equal to a constant g and
  !    f=0, phi=g.x minus its fluid mean solves the discrete problem exactly;
  !    the FPM gradients and densities then have closed forms.
  ! 2. FPM influence field (body_grad=e_x, zero outer data): compatibility
  !    closes, and a distributed solve equals a single-rank host solve of the
  !    same problem (Jacobi: fixed sweeps, to roundoff; PCG: converged).
  ! Both iterations (Jacobi and PCG) are checked.
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use mpi
  use fpm_test_grid
  use fpm_stairstep_mpi_3d
  implicit none
  integer,parameter::gnx=24,gny=20,nz=10
  real(real64),parameter::g(3)=[0.3_real64,-0.7_real64,0.5_real64]
  real(real64)::xf(0:gnx),yf(0:gny),zf(0:nz),gxc(gnx),gyc(gny),zc(nz),gdx(gnx),gdy(gny),dz(nz)
  real(real64),allocatable::phi(:,:,:),phi_ref(:,:,:)
  real(real64)::errors(4),correction,imbalance,residual,check_residual,totals(2),expected(2)
  real(real64)::diff,local,value,ref_totals(2),unused(2)
  character(8),parameter::methods(2)=[character(8)::'jacobi','pcg']
  real(real64)::tolerance,field_tol
  integer::rank,ierr,sweeps,ref_sweeps,ioff,joff,nx,ny,m,max_it
  logical::use_gpu,passed,ok
  character(16)::arg

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  use_gpu=.false.
  if(command_argument_count()>0)then;call get_command_argument(1,arg);use_gpu=trim(arg)=='openacc';endif
  call stretched_axis(gnx,0.0_real64,1.2_real64,0.21_real64,xf,gxc,gdx)
  call stretched_axis(gny,0.0_real64,1.0_real64,-0.17_real64,yf,gyc,gdy)
  call stretched_axis(nz,0.0_real64,0.5_real64,0.13_real64,zf,zc,dz)

  passed=.true.
  do m=1,2
    ! 1. Exact linear solution.
    ! CG reaches a roundoff floor near 1e-11 here; Jacobi goes lower.
    tolerance=merge(1e-12_real64,1e-10_real64,m==1)
    call run_case(MPI_COMM_WORLD,use_gpu,.true.,methods(m),tolerance,400000,phi,ioff,joff,nx,ny, &
         sweeps,residual,check_residual,correction,imbalance,errors,totals,expected)
    ! The final residual is recomputed after the gauge shift, which adds roundoff.
    ok=abs(correction)<1e-13.and.abs(imbalance)<1e-13.and.residual<=2*tolerance.and.check_residual<1e-9
    ok=ok.and.errors(1)<1e-9.and.errors(2)<1e-9.and.errors(3)<1e-9.and.errors(4)==0
    ok=ok.and.maxval(abs(totals-expected))<1e-10
    ! NaN compares false and max() may drop it, so require finite values.
    ok=ok.and.all(ieee_is_finite([residual,check_residual,errors,totals,expected]))
    passed=passed.and.ok
    if(rank==0)then
      write(output_unit,'(a)')'--- method '//trim(methods(m))
      write(output_unit,'(a,2(es12.4,1x))')'exact case compatibility correction/imbalance: ',correction,imbalance
      write(output_unit,'(a,i0,a,2(es12.4,1x))')'exact case iterations: ',sweeps,', residual (solver/independent): ', &
           residual,check_residual
      write(output_unit,'(a,4(es12.4,1x))')'exact case errors phi/grad/vortex/solid (Linf): ',errors
      write(output_unit,'(a,es12.4)')'exact case integral error: ',maxval(abs(totals-expected))
    endif

    ! 2. FPM influence field, distributed versus one host rank: identical fixed
    ! work for Jacobi; converged fields for PCG (its reductions depend on the
    ! decomposition).
    if(m==1)then;tolerance=-1;max_it=600;field_tol=1e-13_real64
    else;tolerance=1e-11_real64;max_it=100000;field_tol=1e-8_real64;endif
    call run_case(MPI_COMM_WORLD,use_gpu,.false.,methods(m),tolerance,max_it,phi,ioff,joff,nx,ny, &
         sweeps,residual,check_residual,correction,imbalance,errors,totals,expected)
    ok=abs(imbalance)<1e-13.and.abs(check_residual-residual)<=1e-8*residual+1e-12
    if(m==1)ok=ok.and.sweeps==600
    if(m==2)ok=ok.and.residual<=2e-11
    call run_case(MPI_COMM_SELF,.false.,.false.,methods(m),tolerance,max_it,phi_ref,ioff,joff,nx,ny, &
         ref_sweeps,value,check_residual,correction,imbalance,errors,ref_totals,unused)
    ! run_case on MPI_COMM_SELF returns the whole grid; compare this rank's block.
    call block_offsets(ioff,joff)
    local=maxval(abs(phi(1:size(phi,1)-2,1:size(phi,2)-2,:)- &
         phi_ref(ioff+1:ioff+size(phi,1)-2,joff+1:joff+size(phi,2)-2,:)))
    call MPI_Allreduce(local,diff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    ok=ok.and.diff<field_tol.and.maxval(abs(totals-ref_totals))<field_tol
    ok=ok.and.all(ieee_is_finite([residual,value,check_residual,diff,totals,ref_totals]))
    ok=ok.and.all(ieee_is_finite(phi(1:size(phi,1)-2,1:size(phi,2)-2,:)))
    if(m==1)ok=ok.and.abs(value-residual)<=1e-8*value
    passed=passed.and.ok
    if(rank==0)then
      write(output_unit,'(a,es12.4)')'influence case compatibility correction: ',correction
      write(output_unit,'(a,2(i0,1x),2(es12.4,1x))')'influence case iterations and residual distributed/one rank: ', &
           sweeps,ref_sweeps,residual,value
      write(output_unit,'(a,es12.4)')'influence case distributed versus one-rank field (Linf): ',diff
      write(output_unit,'(a,2(es22.14,1x))')'influence case vortex/advection integrals: ',totals
      write(output_unit,'(a,es12.4)')'influence case integral difference versus one rank: ',maxval(abs(totals-ref_totals))
    endif
  enddo
  if(.not.passed)then
    if(rank==0)write(error_unit,'(a)')'FAIL: stair-step exactness and decomposition test'
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  endif
  if(rank==0)write(output_unit,'(a,l1,a)')'PASS: stair-step exactness and decomposition test (openacc=',use_gpu,')'
  call MPI_Finalize(ierr)
contains
  pure logical function is_solid(i,j,k)
    integer,intent(in)::i,j,k
    is_solid=(i>=6.and.i<=10.and.j>=5.and.j<=12.and.k>=3.and.k<=7) &
         .or.(i>=15.and.i<=19.and.j>=12.and.j<=16.and.k>=2.and.k<=9) &
         .or.(i>=15.and.i<=17.and.j>=8.and.j<=11.and.k>=2.and.k<=9) &
         .or.(i>=22.and.j<=3)
  end function is_solid

  subroutine block_offsets(ioff,joff)
    integer,intent(out)::ioff,joff
    integer::cart,nx,ny,w,e,s,n
    call cartesian_layout(MPI_COMM_WORLD,gnx,gny,cart,nx,ny,ioff,joff,w,e,s,n)
    call MPI_Comm_free(cart,ierr)
  end subroutine block_offsets

  subroutine run_case(comm,gpu,exact,method,tolerance,max_sweeps,phi,ioff,joff,nx,ny, &
       sweeps,residual,check_residual,correction,imbalance,errors,totals,expected)
    integer,intent(in)::comm,max_sweeps
    logical,intent(in)::gpu,exact
    character(*),intent(in)::method
    real(real64),intent(in)::tolerance
    real(real64),allocatable,intent(inout)::phi(:,:,:)
    integer,intent(out)::ioff,joff,nx,ny,sweeps
    real(real64),intent(out)::residual,check_residual,correction,imbalance,errors(4),totals(2),expected(2)
    integer,allocatable::solid(:,:,:)
    real(real64),allocatable::xc(:),yc(:),dx(:),dy(:),rhs(:,:,:),xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
    real(real64),allocatable::u(:,:,:),v(:,:,:),w(:,:,:),ox(:,:,:),oy(:,:,:),oz(:,:,:),ax(:,:,:),ay(:,:,:),az(:,:,:)
    real(real64),allocatable::gx(:,:,:),gy(:,:,:),gz(:,:,:),vortex(:,:,:),advection(:,:,:)
    real(real64)::body(3),local(2),global(2),mean,bounds(6),ev,cx,cy,cz,vol,lerr(4)
    integer::cart,west,east,south,north,i,j,k,ierr
    call cartesian_layout(comm,gnx,gny,cart,nx,ny,ioff,joff,west,east,south,north)
    allocate(xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny),solid(0:nx+1,0:ny+1,nz),rhs(nx,ny,nz))
    allocate(xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny))
    if(allocated(phi))deallocate(phi)
    allocate(phi(0:nx+1,0:ny+1,nz))
    call local_axis(gnx,gxc,gdx,ioff,nx,xc,dx)
    call local_axis(gny,gyc,gdy,joff,ny,yc,dy)
    solid=0
    do k=1,nz;do j=1,ny;do i=1,nx
      if(is_solid(ioff+i,joff+j,k))solid(i,j,k)=1
    enddo;enddo;enddo
    call exchange_solid_halo_3d(solid,nx,ny,nz,west,east,south,north,cart)
    rhs=0;phi=0
    if(exact)then
      body=g;xlo=g(1);xhi=g(1);ylo=g(2);yhi=g(2);zlo=g(3);zhi=g(3)
    else
      body=[1.0_real64,0.0_real64,0.0_real64];xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0
    endif
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,correction,imbalance)
    call solve_stairstep_poisson_3d(phi,rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,tolerance,max_sweeps,50,sweeps,residual,use_openacc=gpu, &
         method=trim(method))
    check_residual=stairstep_residual_3d(phi,rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart)
    allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz))
    allocate(ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz),gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz))
    allocate(vortex(nx,ny,nz),advection(nx,ny,nz))
    do k=1,nz;do j=1,ny;do i=1,nx
      u(i,j,k)=1+sin(3*xc(i))*cos(2*yc(j));v(i,j,k)=0.5_real64*cos(xc(i)+zc(k));w(i,j,k)=0.2_real64*yc(j)
      ox(i,j,k)=cos(4*zc(k));oy(i,j,k)=-0.3_real64+xc(i)*yc(j);oz(i,j,k)=sin(5*yc(j))
      ax(i,j,k)=0.1_real64*xc(i);ay(i,j,k)=-0.2_real64;az(i,j,k)=zc(k)*yc(j)
    enddo;enddo;enddo
    bounds=[-1.0_real64,2.0_real64,-1.0_real64,2.0_real64,-1.0_real64,1.0_real64]
    call compute_stairstep_fpm_terms_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,solid,body,xc,yc,zc,dx,dy,dz, &
         xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vortex,advection,totals, &
         nx,ny,nz,west,east,south,north,cart,use_openacc=gpu)
    errors=0;expected=0
    if(exact)then
      local=0
      do k=1,nz;do j=1,ny;do i=1,nx
        if(solid(i,j,k)/=0)cycle
        vol=dx(i)*dy(j)*dz(k)
        local=local+[vol*dot_product(g,[xc(i),yc(j),zc(k)]),vol]
      enddo;enddo;enddo
      call MPI_Allreduce(local,global,2,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
      mean=global(1)/global(2)
      lerr=0;local=0
      do k=1,nz;do j=1,ny;do i=1,nx
        if(solid(i,j,k)/=0)then
          lerr(4)=max(lerr(4),abs(phi(i,j,k)),abs(vortex(i,j,k)),abs(gx(i,j,k)))
          cycle
        endif
        lerr(1)=max(lerr(1),abs(phi(i,j,k)-(dot_product(g,[xc(i),yc(j),zc(k)])-mean)))
        lerr(2)=max(lerr(2),abs(gx(i,j,k)-g(1)),abs(gy(i,j,k)-g(2)),abs(gz(i,j,k)-g(3)))
        cx=w(i,j,k)*oy(i,j,k)-v(i,j,k)*oz(i,j,k)
        cy=u(i,j,k)*oz(i,j,k)-w(i,j,k)*ox(i,j,k)
        cz=v(i,j,k)*ox(i,j,k)-u(i,j,k)*oy(i,j,k)
        ev=-(cx*g(1)+cy*g(2)+cz*g(3))
        lerr(3)=max(lerr(3),abs(vortex(i,j,k)-ev))
        vol=dx(i)*dy(j)*dz(k)
        local=local+vol*[ev,-(ax(i,j,k)*g(1)+ay(i,j,k)*g(2)+az(i,j,k)*g(3))]
      enddo;enddo;enddo
      call MPI_Allreduce(lerr,errors,4,MPI_DOUBLE_PRECISION,MPI_MAX,cart,ierr)
      call MPI_Allreduce(local,expected,2,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
    endif
    call MPI_Comm_free(cart,ierr)
  end subroutine run_case
end program test_stairstep_exact_3d
