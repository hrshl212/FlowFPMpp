module fpm_stairstep_mpi_3d
  ! FPM influence-field solve and volume terms around a stair-step (blocked-cell)
  ! immersed body on a stretched Cartesian grid with an x-y domain decomposition.
  !
  ! Solid cells are excluded from the solve and hold phi=0. On every face between
  ! a fluid and a solid cell the axis-aligned face gradient is prescribed by
  ! body_grad(axis). For the influence field of force direction i the body
  ! condition n.grad(phi_i)=n_i gives body_grad=e_i (it is invariant under n->-n).
  ! Outer faces carry prescribed normal derivatives (xlo ... zhi, as d(phi)/dx
  ! etc., not outward-normal values). The pure Neumann problem is made compatible
  ! by one uniform correction of the outer-face data and fixed by removing the
  ! fluid-volume mean. Fields carry one x-y halo layer; z is not decomposed.
  ! use_openacc (default true) selects device execution; false runs the same
  ! code on the host.
  use,intrinsic :: iso_fortran_env,only:real64
  use mpi
  use fpm_halo_mpi_3d
  use fpm_trace
  implicit none
  private
  public::exchange_solid_halo_3d,enforce_stairstep_compatibility_3d,solve_stairstep_poisson_3d
  public::stairstep_residual_3d,compute_stairstep_fpm_terms_3d
contains
  ! Fills the one-layer x-y halo of the solid flag (1 solid, 0 fluid).
  subroutine exchange_solid_halo_3d(solid,nx,ny,nz,west,east,south,north,comm)
    integer,intent(in)::nx,ny,nz,west,east,south,north,comm
    integer,intent(inout)::solid(0:nx+1,0:ny+1,nz)
    type(integer_halo_workspace_3d)::work
    if(west==MPI_PROC_NULL)solid(0,:,:)=0
    if(east==MPI_PROC_NULL)solid(nx+1,:,:)=0
    if(south==MPI_PROC_NULL)solid(:,0,:)=0
    if(north==MPI_PROC_NULL)solid(:,ny+1,:)=0
    call initialize_integer_halo_workspace_3d(work,nx,ny,nz,1)
    call exchange_integer_xy_halo_3d(solid,nx,ny,nz,1,west,east,south,north,comm,work)
    call finalize_integer_halo_workspace_3d(work)
  end subroutine exchange_solid_halo_3d

  ! Shifts the outer-face derivatives by one uniform correction so that the
  ! global outer plus body flux balances the global fluid source integral.
  ! solid must have a current halo.
  subroutine enforce_stairstep_compatibility_3d(rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
       nx,ny,nz,west,east,south,north,comm,correction,imbalance)
    integer,intent(in)::nx,ny,nz,west,east,south,north,comm,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::rhs(nx,ny,nz),body_grad(3),dx(nx),dy(ny),dz(nz)
    real(real64),intent(inout)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(out)::correction,imbalance
    real(real64)::local(3),global(3),source
    integer::i,j,k,ierr
    source=0
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)==0)source=source+rhs(i,j,k)*dx(i)*dy(j)*dz(k)
    enddo;enddo;enddo
    local(1)=source
    call total_flux(local(2),local(3))
    call MPI_Allreduce(local,global,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(global(3)<=0)error stop 'stair-step compatibility needs at least one fluid outer face'
    correction=-(global(2)-global(1))/global(3)
    if(west==MPI_PROC_NULL)xlo=xlo-correction
    if(east==MPI_PROC_NULL)xhi=xhi+correction
    if(south==MPI_PROC_NULL)ylo=ylo-correction
    if(north==MPI_PROC_NULL)yhi=yhi+correction
    zlo=zlo-correction;zhi=zhi+correction
    call total_flux(local(2),local(3))
    call MPI_Allreduce(local(2),global(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    imbalance=global(2)-global(1)
  contains
    ! Net outward flux through fluid outer faces and body faces, and the area
    ! of the fluid outer faces.
    subroutine total_flux(flux,area)
      real(real64),intent(out)::flux,area
      real(real64)::a
      flux=0;area=0
      do k=1,nz;do j=1,ny;do i=1,nx
        if(solid(i,j,k)/=0)cycle
        a=dy(j)*dz(k)
        if(i==1.and.west==MPI_PROC_NULL)then;flux=flux-xlo(j,k)*a;area=area+a
        elseif(solid(i-1,j,k)/=0)then;flux=flux-body_grad(1)*a;endif
        if(i==nx.and.east==MPI_PROC_NULL)then;flux=flux+xhi(j,k)*a;area=area+a
        elseif(solid(i+1,j,k)/=0)then;flux=flux+body_grad(1)*a;endif
        a=dx(i)*dz(k)
        if(j==1.and.south==MPI_PROC_NULL)then;flux=flux-ylo(i,k)*a;area=area+a
        elseif(solid(i,j-1,k)/=0)then;flux=flux-body_grad(2)*a;endif
        if(j==ny.and.north==MPI_PROC_NULL)then;flux=flux+yhi(i,k)*a;area=area+a
        elseif(solid(i,j+1,k)/=0)then;flux=flux+body_grad(2)*a;endif
        a=dx(i)*dy(j)
        if(k==1)then;flux=flux-zlo(i,j)*a;area=area+a
        elseif(solid(i,j,k-1)/=0)then;flux=flux-body_grad(3)*a;endif
        if(k==nz)then;flux=flux+zhi(i,j)*a;area=area+a
        elseif(solid(i,j,k+1)/=0)then;flux=flux+body_grad(3)*a;endif
      enddo;enddo;enddo
    end subroutine total_flux
  end subroutine enforce_stairstep_compatibility_3d

  ! Solves the stair-step pure-Neumann problem over fluid cells, then removes the
  ! global fluid-volume mean. method selects the iteration:
  !   'pcg' (default): conjugate gradients on the volume-weighted (symmetric)
  !     operator with a diagonal preconditioner;
  !   'jacobi': weighted Jacobi, phi <- 0.2 phi + 0.8 (sum(a_nb phi_nb) - f) / sum(a_nb).
  ! Both stop when the max-norm residual of L(phi)-f, checked every
  ! check_interval iterations, is at most tolerance (a negative tolerance runs
  ! exactly max_sweeps iterations). Prescribed outer and body face gradients are
  ! folded once into a source term and solid and outer neighbours are dropped
  ! from the diagonal; solid cells stay at zero, so the kernels read no masks.
  ! solid must have a current halo.
  subroutine solve_stairstep_poisson_3d(phi,rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
       nx,ny,nz,west,east,south,north,comm,tolerance,max_sweeps,check_interval,sweeps,residual,use_openacc,method)
    integer,intent(in)::nx,ny,nz,west,east,south,north,comm,max_sweeps,check_interval
    integer,intent(in)::solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(inout)::phi(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::rhs(nx,ny,nz),body_grad(3)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(in)::xc(0:nx+1),yc(0:ny+1),zc(nz),dx(nx),dy(ny),dz(nz),tolerance
    integer,intent(out)::sweeps
    real(real64),intent(out)::residual
    logical,intent(in),optional::use_openacc
    character(*),intent(in),optional::method
    type(halo_workspace_3d)::work
    real(real64),allocatable::next(:,:,:),src(:,:,:),diag(:,:,:),rdiag(:,:,:),r(:,:,:),q(:,:,:)
    real(real64),allocatable::cw(:),ce(:),cs(:),cn(:),cb(:),ct(:)
    real(real64)::local(2),global(2),mean_value
    integer::i,j,k,ierr,isolated
    logical::accelerate,use_pcg
    accelerate=.true.;if(present(use_openacc))accelerate=use_openacc
    use_pcg=.true.
    if(present(method))then
      if(method/='pcg'.and.method/='jacobi')error stop 'stair-step method must be pcg or jacobi'
      use_pcg=method=='pcg'
    endif
    if(nx<1.or.ny<1.or.nz<2)error stop 'stair-step solve needs nz>=2 and nonempty ranks'
    call trace_push('stair-step setup')
    ! next is the second Jacobi buffer or the PCG search direction.
    allocate(next(0:nx+1,0:ny+1,nz),src(nx,ny,nz),diag(nx,ny,nz),rdiag(nx,ny,nz))
    if(use_pcg)then;allocate(r(nx,ny,nz),q(nx,ny,nz));else;allocate(r(1,1,1),q(1,1,1));endif
    allocate(cw(nx),ce(nx),cs(ny),cn(ny),cb(nz),ct(nz))
    call axis_coefficients(xc,dx,nx,west,east,cw,ce)
    call axis_coefficients(yc,dy,ny,south,north,cs,cn)
    do k=1,nz
      cb(k)=0;ct(k)=0
      if(k>1)cb(k)=1/(dz(k)*(zc(k)-zc(k-1)))
      if(k<nz)ct(k)=1/(dz(k)*(zc(k+1)-zc(k)))
    enddo
    call initialize_halo_workspace_3d(work,nx,ny,nz,1)
    sweeps=0;residual=huge(1.0_real64)
    ! Only phi moves between host and device (and is skipped if the caller
    ! already holds it on the device); the operator arrays are built in place.
    !$acc data copy(phi) copyin(rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,cw,ce,cs,cn,cb,ct,dx,dy,dz) &
    !$acc& create(next,src,diag,rdiag,r,q) if(accelerate)
    call build_operator(rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,cw,ce,cs,cn,cb,ct,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,use_pcg,src,diag,rdiag,isolated,accelerate)
    if(isolated>0)error stop 'stair-step fluid cell has no fluid neighbour'
    ! Solid cells and physical-boundary halos hold zero in phi and next; they
    ! carry zero weight but must be finite.
    !$acc parallel loop collapse(3) default(present) if(accelerate)
    do k=1,nz;do j=0,ny+1;do i=0,nx+1
      if(solid(i,j,k)/=0.or.(i==0.and.west==MPI_PROC_NULL).or.(i==nx+1.and.east==MPI_PROC_NULL) &
           .or.(j==0.and.south==MPI_PROC_NULL).or.(j==ny+1.and.north==MPI_PROC_NULL))phi(i,j,k)=0
      next(i,j,k)=phi(i,j,k)
    enddo;enddo;enddo
    !$acc end parallel loop
    call trace_pop()
    call trace_push('stair-step iterations')
    if(use_pcg)then
      call pcg_iterate()
    else
      call jacobi_iterate()
    endif
    call trace_pop()
    call trace_push('stair-step gauge')
    local=0
    !$acc parallel loop collapse(3) present(phi,dx,dy,dz,solid) reduction(+:local) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)==0)then
        local(1)=local(1)+phi(i,j,k)*dx(i)*dy(j)*dz(k)
        local(2)=local(2)+dx(i)*dy(j)*dz(k)
      endif
    enddo;enddo;enddo
    !$acc end parallel loop
    call MPI_Allreduce(local,global,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    mean_value=global(1)/global(2)
    !$acc parallel loop collapse(3) present(phi,solid) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)==0)phi(i,j,k)=phi(i,j,k)-mean_value
    enddo;enddo;enddo
    !$acc end parallel loop
    call exchange(phi)
    residual=operator_residual(phi,src,diag,cw,ce,cs,cn,cb,ct,nx,ny,nz,comm,accelerate)
    call trace_pop()
    !$acc end data
    call finalize_halo_workspace_3d(work)
    deallocate(next,src,diag,rdiag,r,q,cw,ce,cs,cn,cb,ct)
  contains
    subroutine exchange(field)
      real(real64),intent(inout)::field(0:nx+1,0:ny+1,nz)
      call exchange_xy_halo_3d(field,nx,ny,nz,1,west,east,south,north,comm,work, &
           host_only=.not.accelerate,batched_device=accelerate)
    end subroutine exchange

    ! Ping-pong between phi and next; an odd step count ends with one copy.
    subroutine jacobi_iterate()
      integer::step,nsteps,i,j,k,n1,n2,n3
      n1=nx;n2=ny;n3=nz
      do while(sweeps<max_sweeps.and.residual>tolerance)
        nsteps=min(max(1,check_interval),max_sweeps-sweeps)
        do step=1,nsteps
          if(mod(step,2)==1)then
            call exchange(phi)
            call jacobi_sweep(phi,next,src,rdiag,cw,ce,cs,cn,cb,ct,nx,ny,nz,accelerate)
          else
            call exchange(next)
            call jacobi_sweep(next,phi,src,rdiag,cw,ce,cs,cn,cb,ct,nx,ny,nz,accelerate)
          endif
          sweeps=sweeps+1
        enddo
        if(mod(nsteps,2)==1)then
          !$acc parallel loop collapse(3) present(phi,next) if(accelerate)
          do k=1,n3;do j=1,n2;do i=1,n1
            phi(i,j,k)=next(i,j,k)
          enddo;enddo;enddo
          !$acc end parallel loop
        endif
        call exchange(phi)
        residual=operator_residual(phi,src,diag,cw,ce,cs,cn,cb,ct,nx,ny,nz,comm,accelerate)
      enddo
    end subroutine jacobi_iterate

    ! Preconditioned CG on A phi = b with A = V (diag - offdiag) and b = V src,
    ! where V is the cell volume. A is symmetric positive semidefinite on the
    ! fluid cells with the constants as null space; compatibility makes b
    ! orthogonal to it, and every true residual is projected once more to remove
    ! roundoff. rdiag holds 1/(V diag) (the preconditioner) and next the search
    ! direction p. At each convergence check the recursively updated residual is
    ! replaced by the true one (it drifts near roundoff), and a breakdown
    ! restarts from the steepest-descent direction.
    subroutine pcg_iterate()
      real(real64)::rz,rz_new,pq,alpha,beta,vol,res
      integer::it,i,j,k,ierr,n1,n2,n3,restarts
      ! The residual is passed back explicitly: NVHPC OpenACC builds did not see
      ! an update of the host-associated residual made in a sibling procedure.
      n1=nx;n2=ny;n3=nz
      call true_residual(.true.,rz,res)
      it=0;restarts=0
      do while(sweeps<max_sweeps.and.res>tolerance)
        call exchange(next)
        pq=0
        !$acc parallel loop collapse(3) present(next,q,diag,cw,ce,cs,cn,cb,ct,dx,dy,dz) private(vol) &
        !$acc& reduction(+:pq) if(accelerate)
        do k=1,n3;do j=1,n2;do i=1,n1
          vol=dx(i)*dy(j)*dz(k)
          q(i,j,k)=merge(vol*(diag(i,j,k)*next(i,j,k)-cw(i)*next(i-1,j,k)-ce(i)*next(i+1,j,k) &
               -cs(j)*next(i,j-1,k)-cn(j)*next(i,j+1,k)-cb(k)*next(i,j,max(k-1,1)) &
               -ct(k)*next(i,j,min(k+1,n3))),0.0_real64,diag(i,j,k)>0)
          pq=pq+next(i,j,k)*q(i,j,k)
        enddo;enddo;enddo
        !$acc end parallel loop
        call MPI_Allreduce(MPI_IN_PLACE,pq,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        if(pq<=0.or.rz<=0)then
          ! Breakdown (converged to roundoff or lost conjugacy): restart once.
          restarts=restarts+1
          if(restarts>1)exit
          call true_residual(.true.,rz,res)
          cycle
        endif
        restarts=0
        alpha=rz/pq
        rz_new=0
        !$acc parallel loop collapse(3) present(phi,next,r,q,rdiag) reduction(+:rz_new) if(accelerate)
        do k=1,n3;do j=1,n2;do i=1,n1
          phi(i,j,k)=phi(i,j,k)+alpha*next(i,j,k)
          r(i,j,k)=r(i,j,k)-alpha*q(i,j,k)
          rz_new=rz_new+rdiag(i,j,k)*r(i,j,k)*r(i,j,k)
        enddo;enddo;enddo
        !$acc end parallel loop
        call MPI_Allreduce(MPI_IN_PLACE,rz_new,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        beta=rz_new/rz;rz=rz_new
        !$acc parallel loop collapse(3) present(next,r,rdiag) if(accelerate)
        do k=1,n3;do j=1,n2;do i=1,n1
          next(i,j,k)=rdiag(i,j,k)*r(i,j,k)+beta*next(i,j,k)
        enddo;enddo;enddo
        !$acc end parallel loop
        sweeps=sweeps+1;it=it+1
        if(mod(it,max(1,check_interval))==0.or.sweeps==max_sweeps)call true_residual(.false.,rz,res)
      enddo
      residual=res
    end subroutine pcg_iterate

    ! Sets r to the projected true residual b - A phi, rz to r.M^-1 r, and res
    ! to the max-norm residual of L(phi)-f; reset also sets p = M^-1 r.
    subroutine true_residual(reset,rz,res)
      logical,intent(in)::reset
      real(real64),intent(out)::rz,res
      real(real64)::red(2),gred(2),vol,shift
      integer::i,j,k,ierr,n1,n2,n3
      n1=nx;n2=ny;n3=nz
      call exchange(phi)
      res=operator_residual(phi,src,diag,cw,ce,cs,cn,cb,ct,nx,ny,nz,comm,accelerate)
      red=0
      !$acc parallel loop collapse(3) present(phi,src,diag,cw,ce,cs,cn,cb,ct,dx,dy,dz,r) private(vol) &
      !$acc& reduction(+:red) if(accelerate)
      do k=1,n3;do j=1,n2;do i=1,n1
        r(i,j,k)=0
        if(diag(i,j,k)>0)then
          vol=dx(i)*dy(j)*dz(k)
          r(i,j,k)=vol*(src(i,j,k)+cw(i)*phi(i-1,j,k)+ce(i)*phi(i+1,j,k)+cs(j)*phi(i,j-1,k) &
               +cn(j)*phi(i,j+1,k)+cb(k)*phi(i,j,max(k-1,1))+ct(k)*phi(i,j,min(k+1,n3))-diag(i,j,k)*phi(i,j,k))
          red(1)=red(1)+r(i,j,k);red(2)=red(2)+1
        endif
      enddo;enddo;enddo
      !$acc end parallel loop
      call MPI_Allreduce(red,gred,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      shift=gred(1)/gred(2)
      rz=0
      !$acc parallel loop collapse(3) present(r,rdiag,diag,next) reduction(+:rz) if(accelerate)
      do k=1,n3;do j=1,n2;do i=1,n1
        if(diag(i,j,k)>0)r(i,j,k)=r(i,j,k)-shift
        if(reset)next(i,j,k)=rdiag(i,j,k)*r(i,j,k)
        rz=rz+rdiag(i,j,k)*r(i,j,k)*r(i,j,k)
      enddo;enddo;enddo
      !$acc end parallel loop
      call MPI_Allreduce(MPI_IN_PLACE,rz,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    end subroutine true_residual
  end subroutine solve_stairstep_poisson_3d

  ! Source term, diagonal, and inverse-diagonal factor of the stair-step
  ! operator (all zero in solid cells), built where the arrays live. rdiag is
  ! 0.8/diag for Jacobi and 1/(V diag) for PCG. isolated counts fluid cells
  ! without any fluid neighbour.
  subroutine build_operator(rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,cw,ce,cs,cn,cb,ct,dx,dy,dz, &
       nx,ny,nz,west,east,south,north,for_pcg,src,diag,rdiag,isolated,accelerate)
    integer,intent(in)::nx,ny,nz,west,east,south,north,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::rhs(nx,ny,nz),body_grad(3)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(in)::cw(nx),ce(nx),cs(ny),cn(ny),cb(nz),ct(nz),dx(nx),dy(ny),dz(nz)
    logical,intent(in)::for_pcg,accelerate
    real(real64),intent(out)::src(nx,ny,nz),diag(nx,ny,nz),rdiag(nx,ny,nz)
    integer,intent(out)::isolated
    real(real64)::s,d
    integer::i,j,k
    isolated=0
    !$acc parallel loop collapse(3) default(present) private(s,d) reduction(+:isolated) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      s=0;d=0
      if(solid(i,j,k)==0)then
        s=-rhs(i,j,k)
        if(i==1.and.west==MPI_PROC_NULL)then;s=s-xlo(j,k)/dx(i)
        elseif(solid(i-1,j,k)/=0)then;s=s-body_grad(1)/dx(i)
        else;d=d+cw(i);endif
        if(i==nx.and.east==MPI_PROC_NULL)then;s=s+xhi(j,k)/dx(i)
        elseif(solid(i+1,j,k)/=0)then;s=s+body_grad(1)/dx(i)
        else;d=d+ce(i);endif
        if(j==1.and.south==MPI_PROC_NULL)then;s=s-ylo(i,k)/dy(j)
        elseif(solid(i,j-1,k)/=0)then;s=s-body_grad(2)/dy(j)
        else;d=d+cs(j);endif
        if(j==ny.and.north==MPI_PROC_NULL)then;s=s+yhi(i,k)/dy(j)
        elseif(solid(i,j+1,k)/=0)then;s=s+body_grad(2)/dy(j)
        else;d=d+cn(j);endif
        if(k==1)then;s=s-zlo(i,j)/dz(k)
        elseif(solid(i,j,k-1)/=0)then;s=s-body_grad(3)/dz(k)
        else;d=d+cb(k);endif
        if(k==nz)then;s=s+zhi(i,j)/dz(k)
        elseif(solid(i,j,k+1)/=0)then;s=s+body_grad(3)/dz(k)
        else;d=d+ct(k);endif
        if(d<=0)isolated=isolated+1
      endif
      src(i,j,k)=s;diag(i,j,k)=d
      if(d<=0)then;rdiag(i,j,k)=0
      elseif(for_pcg)then;rdiag(i,j,k)=1/(dx(i)*dy(j)*dz(k)*d)
      else;rdiag(i,j,k)=0.8_real64/d;endif
    enddo;enddo;enddo
    !$acc end parallel loop
  end subroutine build_operator

  ! One weighted-Jacobi sweep from a (with current halos) into b. Solid cells
  ! have src=rdiag=0 and a=0, so they stay zero; solid and outer neighbours
  ! contribute nothing because they hold zero.
  subroutine jacobi_sweep(a,b,src,rdiag,cw,ce,cs,cn,cb,ct,nx,ny,nz,accelerate)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(in)::a(0:nx+1,0:ny+1,nz),src(nx,ny,nz),rdiag(nx,ny,nz)
    real(real64),intent(in)::cw(nx),ce(nx),cs(ny),cn(ny),cb(nz),ct(nz)
    real(real64),intent(inout)::b(0:nx+1,0:ny+1,nz)
    logical,intent(in)::accelerate
    integer::i,j,k
    !$acc parallel loop collapse(3) present(a,b,src,rdiag,cw,ce,cs,cn,cb,ct) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      ! Damping removes the checkerboard two-cycle of undamped Jacobi.
      b(i,j,k)=0.2_real64*a(i,j,k)+rdiag(i,j,k)* &
           (src(i,j,k)+cw(i)*a(i-1,j,k)+ce(i)*a(i+1,j,k)+cs(j)*a(i,j-1,k)+cn(j)*a(i,j+1,k) &
           +cb(k)*a(i,j,max(k-1,1))+ct(k)*a(i,j,min(k+1,nz)))
    enddo;enddo;enddo
    !$acc end parallel loop
  end subroutine jacobi_sweep

  ! Independent host evaluation of the global max-norm residual L(phi)-f over
  ! fluid cells, from face gradients. phi and solid halos must be current.
  function stairstep_residual_3d(phi,rhs,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
       nx,ny,nz,west,east,south,north,comm) result(value)
    integer,intent(in)::nx,ny,nz,west,east,south,north,comm,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,nz),rhs(nx,ny,nz),body_grad(3)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(in)::xc(0:nx+1),yc(0:ny+1),zc(nz),dx(nx),dy(ny),dz(nz)
    real(real64)::value,local,g(6)
    integer::i,j,k,ierr
    local=0
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)/=0)cycle
      call face_gradients(phi,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,nx,ny,nz, &
           west,east,south,north,i,j,k,g)
      local=max(local,abs((g(2)-g(1))/dx(i)+(g(4)-g(3))/dy(j)+(g(6)-g(5))/dz(k)-rhs(i,j,k)))
    enddo;enddo;enddo
    call MPI_Allreduce(local,value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  end function stairstep_residual_3d

  ! FPM volume densities in fluid cells from grad(phi) (average of the two face
  ! gradients per axis, with prescribed outer and body face gradients):
  ! vortex = -(omega x u).grad(phi), advection = -a.grad(phi), and their global
  ! integrals over fluid cells with centres strictly inside bounds. Solid cells
  ! get zero. phi and solid halos must be current.
  subroutine compute_stairstep_fpm_terms_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,solid,body_grad, &
       xc,yc,zc,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vortex,advection,totals, &
       nx,ny,nz,west,east,south,north,comm,use_openacc)
    integer,intent(in)::nx,ny,nz,west,east,south,north,comm,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
    real(real64),intent(in)::ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz),ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
    real(real64),intent(in)::body_grad(3),xc(0:nx+1),yc(0:ny+1),zc(nz),dx(nx),dy(ny),dz(nz)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny),bounds(6)
    real(real64),intent(out)::gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz),vortex(nx,ny,nz),advection(nx,ny,nz)
    real(real64),intent(out)::totals(2)
    logical,intent(in),optional::use_openacc
    real(real64)::g(6),cx,cy,cz,mask,volume,vortex_total,advection_total,local(2)
    integer::i,j,k,ierr
    logical::accelerate
    accelerate=.true.;if(present(use_openacc))accelerate=use_openacc
    call trace_push('stair-step FPM terms')
    vortex_total=0;advection_total=0
    !$acc data copyin(phi,u,v,w,ox,oy,oz,ax,ay,az,solid,body_grad,xc,yc,zc,dx,dy,dz) &
    !$acc& copyin(xlo,xhi,ylo,yhi,zlo,zhi,bounds) copyout(gx,gy,gz,vortex,advection) if(accelerate)
    !$acc parallel loop collapse(3) default(present) private(g,cx,cy,cz,mask,volume) &
    !$acc& reduction(+:vortex_total,advection_total) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      call face_gradients(phi,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,nx,ny,nz, &
           west,east,south,north,i,j,k,g)
      mask=merge(1.0_real64,0.0_real64,solid(i,j,k)==0)
      gx(i,j,k)=mask*0.5_real64*(g(1)+g(2))
      gy(i,j,k)=mask*0.5_real64*(g(3)+g(4))
      gz(i,j,k)=mask*0.5_real64*(g(5)+g(6))
      cx=w(i,j,k)*oy(i,j,k)-v(i,j,k)*oz(i,j,k)
      cy=u(i,j,k)*oz(i,j,k)-w(i,j,k)*ox(i,j,k)
      cz=v(i,j,k)*ox(i,j,k)-u(i,j,k)*oy(i,j,k)
      vortex(i,j,k)=-(cx*gx(i,j,k)+cy*gy(i,j,k)+cz*gz(i,j,k))
      advection(i,j,k)=-(ax(i,j,k)*gx(i,j,k)+ay(i,j,k)*gy(i,j,k)+az(i,j,k)*gz(i,j,k))
      volume=merge(dx(i)*dy(j)*dz(k),0.0_real64,xc(i)>bounds(1).and.xc(i)<bounds(2).and. &
           yc(j)>bounds(3).and.yc(j)<bounds(4).and.zc(k)>bounds(5).and.zc(k)<bounds(6))
      vortex_total=vortex_total+volume*vortex(i,j,k)
      advection_total=advection_total+volume*advection(i,j,k)
    enddo;enddo;enddo
    !$acc end parallel loop
    !$acc end data
    local=[vortex_total,advection_total]
    call MPI_Allreduce(local,totals,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call trace_pop()
  end subroutine compute_stairstep_fpm_terms_3d

  ! West, east, south, north, bottom, top face gradients of fluid cell (i,j,k).
  subroutine face_gradients(phi,solid,body_grad,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,nx,ny,nz, &
       west,east,south,north,i,j,k,g)
    !$acc routine seq
    integer,intent(in)::nx,ny,nz,west,east,south,north,i,j,k,solid(0:nx+1,0:ny+1,nz)
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,nz),body_grad(3),xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz)
    real(real64),intent(in)::zlo(nx,ny),zhi(nx,ny),xc(0:nx+1),yc(0:ny+1),zc(nz)
    real(real64),intent(out)::g(6)
    if(i==1.and.west==MPI_PROC_NULL)then;g(1)=xlo(j,k)
    elseif(solid(i-1,j,k)/=0)then;g(1)=body_grad(1)
    else;g(1)=(phi(i,j,k)-phi(i-1,j,k))/(xc(i)-xc(i-1));endif
    if(i==nx.and.east==MPI_PROC_NULL)then;g(2)=xhi(j,k)
    elseif(solid(i+1,j,k)/=0)then;g(2)=body_grad(1)
    else;g(2)=(phi(i+1,j,k)-phi(i,j,k))/(xc(i+1)-xc(i));endif
    if(j==1.and.south==MPI_PROC_NULL)then;g(3)=ylo(i,k)
    elseif(solid(i,j-1,k)/=0)then;g(3)=body_grad(2)
    else;g(3)=(phi(i,j,k)-phi(i,j-1,k))/(yc(j)-yc(j-1));endif
    if(j==ny.and.north==MPI_PROC_NULL)then;g(4)=yhi(i,k)
    elseif(solid(i,j+1,k)/=0)then;g(4)=body_grad(2)
    else;g(4)=(phi(i,j+1,k)-phi(i,j,k))/(yc(j+1)-yc(j));endif
    if(k==1)then;g(5)=zlo(i,j)
    elseif(solid(i,j,k-1)/=0)then;g(5)=body_grad(3)
    else;g(5)=(phi(i,j,k)-phi(i,j,k-1))/(zc(k)-zc(k-1));endif
    if(k==nz)then;g(6)=zhi(i,j)
    elseif(solid(i,j,k+1)/=0)then;g(6)=body_grad(3)
    else;g(6)=(phi(i,j,k+1)-phi(i,j,k))/(zc(k+1)-zc(k));endif
  end subroutine face_gradients

  subroutine axis_coefficients(center,width,n,lower,upper,cl,cu)
    integer,intent(in)::n,lower,upper
    real(real64),intent(in)::center(0:n+1),width(n)
    real(real64),intent(out)::cl(n),cu(n)
    integer::m
    do m=1,n
      cl(m)=1/(width(m)*(center(m)-center(m-1)))
      cu(m)=1/(width(m)*(center(m+1)-center(m)))
    enddo
    if(lower==MPI_PROC_NULL)cl(1)=0
    if(upper==MPI_PROC_NULL)cu(n)=0
  end subroutine axis_coefficients

  ! Max-norm of src + sum(a_nb phi_nb) - diag phi over fluid cells (diag>0),
  ! which equals L(phi)-f because src holds -f plus the prescribed face terms.
  function operator_residual(phi,src,diag,cw,ce,cs,cn,cb,ct,nx,ny,nz,comm,accelerate) result(value)
    integer,intent(in)::nx,ny,nz,comm
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,nz),src(nx,ny,nz),diag(nx,ny,nz)
    real(real64),intent(in)::cw(nx),ce(nx),cs(ny),cn(ny),cb(nz),ct(nz)
    logical,intent(in)::accelerate
    real(real64)::value,local
    integer::i,j,k,ierr
    local=0
    !$acc parallel loop collapse(3) present(phi,src,diag,cw,ce,cs,cn,cb,ct) reduction(max:local) if(accelerate)
    do k=1,nz;do j=1,ny;do i=1,nx
      if(diag(i,j,k)>0)local=max(local,abs(src(i,j,k)+cw(i)*phi(i-1,j,k)+ce(i)*phi(i+1,j,k) &
           +cs(j)*phi(i,j-1,k)+cn(j)*phi(i,j+1,k)+cb(k)*phi(i,j,max(k-1,1))+ct(k)*phi(i,j,min(k+1,nz)) &
           -diag(i,j,k)*phi(i,j,k)))
    enddo;enddo;enddo
    !$acc end parallel loop
    call MPI_Allreduce(local,value,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  end function operator_residual
end module fpm_stairstep_mpi_3d
