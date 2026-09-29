program benchmark_stairstep_3d
  ! Public fixed-work benchmark: the x influence field of a stair-step sphere on
  ! a mildly stretched grid (compatibility + Jacobi sweeps + gauge), then the
  ! FPM volume terms from synthetic velocity and vorticity fields.
  ! Usage: benchmark_stairstep_3d [mode [gnx gny nz [iterations [method [tolerance]]]]]
  !   method jacobi (default, fixed-work bandwidth study) or pcg
  !   tolerance > 0: time to solution (stop at that max-norm residual, checked
  !   every 100 iterations, with iterations as the cap); default: fixed work
  !   mode benchmark: CPU and GPU runs with a parity check (default)
  !   mode cpu:       CPU only (any rank count whose layout divides the grid)
  !   mode gpu:       GPU only (for profiling)
  ! Times are medians of three runs after one warm-up and exclude input
  ! generation. On the GPU one data region spans the frame: the flow fields
  ! and phi are uploaded once, phi stays on the device between the solve and
  ! the force terms, and only the two force densities and phi are downloaded.
  ! The solve time includes the uploads and the terms time the downloads.
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use mpi
  use fpm_test_grid
  use fpm_stairstep_mpi_3d
  use fpm_trace
  implicit none
  integer,parameter::repeats=3
  real(real64),allocatable::gxf(:),gxc(:),gdx(:),gyf(:),gyc(:),gdy(:),zf(:),zc(:),dz(:)
  real(real64),allocatable::xc(:),yc(:),dx(:),dy(:),phi(:,:,:),phi_cpu(:,:,:),rhs(:,:,:)
  real(real64),allocatable::xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  real(real64),allocatable::u(:,:,:),v(:,:,:),w(:,:,:),ox(:,:,:),oy(:,:,:),oz(:,:,:),ax(:,:,:),ay(:,:,:),az(:,:,:)
  real(real64),allocatable::gx(:,:,:),gy(:,:,:),gz(:,:,:),vortex(:,:,:),advection(:,:,:)
  integer,allocatable::solid(:,:,:)
  real(real64)::body(3),bounds(6),lz,radius,center(3),t_solve(repeats,2),t_terms(repeats,2)
  real(real64)::totals(2),totals_cpu(2),residual,local(3),global(3),linf,rel_l2,setup_start,setup_seconds
  real(real64)::tolerance
  integer::gnx,gny,nz,nsweep,rank,nproc,ierr,cart,nx,ny,ioff,joff,west,east,south,north,i,j,k,r,sweeps
  integer(8)::nsolid,nsolid_local
  character(16)::mode,arg,method
  mode='benchmark';gnx=256;gny=256;nz=64;nsweep=250;method='jacobi';tolerance=-1
  if(command_argument_count()>0)call get_command_argument(1,mode)
  if(all(trim(mode)/=[character(16)::'benchmark','cpu','gpu']))error stop 'mode must be benchmark, cpu, or gpu'
  if(command_argument_count()>=4)then
    call get_command_argument(2,arg);read(arg,*)gnx
    call get_command_argument(3,arg);read(arg,*)gny
    call get_command_argument(4,arg);read(arg,*)nz
  endif
  if(command_argument_count()>=5)then;call get_command_argument(5,arg);read(arg,*)nsweep;endif
  if(command_argument_count()>=6)call get_command_argument(6,method)
  if(command_argument_count()>=7)then;call get_command_argument(7,arg);read(arg,*)tolerance;endif
  call MPI_Init(ierr)
  setup_start=MPI_Wtime()
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call cartesian_layout(MPI_COMM_WORLD,gnx,gny,cart,nx,ny,ioff,joff,west,east,south,north)
  lz=real(nz,real64)/gnx
  allocate(gxf(0:gnx),gxc(gnx),gdx(gnx),gyf(0:gny),gyc(gny),gdy(gny),zf(0:nz),zc(nz),dz(nz))
  call stretched_axis(gnx,0.0_real64,1.0_real64,0.2_real64,gxf,gxc,gdx)
  call stretched_axis(gny,0.0_real64,real(gny,real64)/gnx,-0.2_real64,gyf,gyc,gdy)
  call stretched_axis(nz,0.0_real64,lz,0.1_real64,zf,zc,dz)
  allocate(xc(0:nx+1),yc(0:ny+1),dx(nx),dy(ny),solid(0:nx+1,0:ny+1,nz))
  allocate(phi(0:nx+1,0:ny+1,nz),phi_cpu(0:nx+1,0:ny+1,nz),rhs(nx,ny,nz))
  allocate(xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny))
  allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz))
  allocate(ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz),gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz))
  allocate(vortex(nx,ny,nz),advection(nx,ny,nz))
  call local_axis(gnx,gxc,gdx,ioff,nx,xc,dx)
  call local_axis(gny,gyc,gdy,joff,ny,yc,dy)
  center=[0.5_real64,0.5_real64*gny/gnx,0.5_real64*lz];radius=0.4_real64*lz
  solid=0
  do k=1,nz;do j=1,ny;do i=1,nx
    if(sum(([xc(i),yc(j),zc(k)]-center)**2)<radius**2)solid(i,j,k)=1
    u(i,j,k)=1+0.1_real64*sin(6*xc(i))*cos(4*yc(j));v(i,j,k)=0.2_real64*cos(3*xc(i)+2*zc(k))
    w(i,j,k)=0.1_real64*sin(5*yc(j));ox(i,j,k)=cos(4*zc(k)+yc(j));oy(i,j,k)=sin(3*xc(i))
    oz(i,j,k)=cos(2*xc(i))*sin(5*yc(j));ax(i,j,k)=0.1_real64*u(i,j,k);ay(i,j,k)=-0.2_real64*v(i,j,k)
    az(i,j,k)=0.3_real64*w(i,j,k)
  enddo;enddo;enddo
  call exchange_solid_halo_3d(solid,nx,ny,nz,west,east,south,north,cart)
  nsolid_local=count(solid(1:nx,1:ny,:)/=0)
  call MPI_Allreduce(nsolid_local,nsolid,1,MPI_INTEGER8,MPI_SUM,cart,ierr)
  body=[1.0_real64,0.0_real64,0.0_real64];rhs=0
  bounds=[-1.0_real64,2.0_real64,-1.0_real64,2.0_real64,-1.0_real64,2.0_real64]
  call MPI_Barrier(cart,ierr)
  setup_seconds=MPI_Wtime()-setup_start
  if(rank==0)then
    write(output_unit,'(a,3(i0,1x),a,i0,a)')'grid nx ny nz: ',gnx,gny,nz,'(',int(gnx,8)*gny*nz,' cells)'
    write(output_unit,'(a,i0,a,i0)')'MPI ranks: ',nproc,', solid cells: ',nsolid
    if(tolerance>0)then
      write(output_unit,'(a,a,a,es10.3,a,i0)')'method ',trim(method),', tolerance ',tolerance,', iteration cap ',nsweep
    else
      write(output_unit,'(a,a,a,i0)')'method ',trim(method),', fixed iterations: ',nsweep
    endif
    write(output_unit,'(a,f10.4)')'input setup seconds: ',setup_seconds
  endif
  t_solve=0;t_terms=0
  if(trim(mode)/='gpu')then
    call run_case(.false.,t_solve(1,1),t_terms(1,1))
    do r=1,repeats;call run_case(.false.,t_solve(r,1),t_terms(r,1));enddo
    phi_cpu=phi;totals_cpu=totals
  endif
  if(trim(mode)/='cpu')then
    call run_case(.true.,t_solve(1,2),t_terms(1,2),'GPU warm-up')
    do r=1,repeats;call run_case(.true.,t_solve(r,2),t_terms(r,2));enddo
  endif
  local(1)=sum(phi(1:nx,1:ny,:));local(2)=sum(abs(phi(1:nx,1:ny,:)))
  call MPI_Allreduce(local(1:2),global(1:2),2,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
  if(rank==0)then
    write(output_unit,'(a,i0,a,es14.6)')'last-run iterations: ',sweeps,', final residual: ',residual
    write(output_unit,'(a,2(es25.17,1x))')'last-run phi checksum (sum, sum abs): ',global(1:2)
    write(output_unit,'(a,2(es22.14,1x))')'vortex/advection integrals: ',totals
    if(trim(mode)/='gpu')write(output_unit,'(a,2(f12.6,1x))')'CPU solve/terms medians (s): ', &
         median3(t_solve(:,1)),median3(t_terms(:,1))
    if(trim(mode)/='cpu')write(output_unit,'(a,2(f12.6,1x))')'GPU solve/terms medians (s): ', &
         median3(t_solve(:,2)),median3(t_terms(:,2))
    if(trim(mode)/='gpu')write(output_unit,'(a,f12.6)')'CPU frame (solve+terms) median (s): ', &
         median3(t_solve(:,1)+t_terms(:,1))
    if(trim(mode)/='cpu')write(output_unit,'(a,f12.6)')'GPU frame (solve+terms) median (s): ', &
         median3(t_solve(:,2)+t_terms(:,2))
  endif
  if(trim(mode)=='benchmark')then
    local=0
    do k=1,nz;do j=1,ny;do i=1,nx
      if(solid(i,j,k)/=0)cycle
      local(1)=max(local(1),abs(phi(i,j,k)-phi_cpu(i,j,k)))
      local(2)=local(2)+(phi(i,j,k)-phi_cpu(i,j,k))**2;local(3)=local(3)+phi_cpu(i,j,k)**2
    enddo;enddo;enddo
    call MPI_Allreduce(local(1),global(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,cart,ierr)
    call MPI_Allreduce(local(2:3),global(2:3),2,MPI_DOUBLE_PRECISION,MPI_SUM,cart,ierr)
    linf=global(1);rel_l2=sqrt(global(2)/max(global(3),tiny(1.0_real64)))
    if(rank==0)then
      write(output_unit,'(a,3(es12.4,1x))')'CPU/GPU fluid phi Linf, relative L2, integral difference: ', &
           linf,rel_l2,maxval(abs(totals-totals_cpu))
      write(output_unit,'(a,2(f8.2,1x))')'GPU speedup solve/terms: ', &
           median3(t_solve(:,1))/median3(t_solve(:,2)),median3(t_terms(:,1))/median3(t_terms(:,2))
    endif
    ! Fixed work is the same arithmetic on both; converged runs agree to the tolerance.
    if(linf>merge(1e-10_real64,1e-6_real64,tolerance<=0).or.rel_l2>merge(1e-10_real64,1e-6_real64,tolerance<=0) &
         .or.maxval(abs(totals-totals_cpu))>merge(1e-10_real64,1e-6_real64,tolerance<=0))then
      if(rank==0)write(error_unit,'(a)')'FAIL: CPU/GPU parity'
      call MPI_Abort(cart,3,ierr)
    endif
    if(rank==0)write(output_unit,'(a)')'PASS: stair-step benchmark CPU/GPU parity'
  endif
  call MPI_Comm_free(cart,ierr)
  call MPI_Finalize(ierr)
contains
  subroutine run_case(gpu,solve_seconds,terms_seconds,label)
    logical,intent(in)::gpu
    real(real64),intent(out)::solve_seconds,terms_seconds
    character(*),intent(in),optional::label
    real(real64)::t0,t1,t2,local_t(2),global_t(2),correction,imbalance
    phi=0;xlo=0;xhi=0;ylo=0;yhi=0;zlo=0;zhi=0
    call MPI_Barrier(cart,ierr)
    if(present(label))then;call trace_push(label);else;call trace_push(merge('GPU solve','CPU solve',gpu));endif
    t0=MPI_Wtime()
    !$acc data copyin(u,v,w,ox,oy,oz,ax,ay,az,solid,rhs) copy(phi) create(gx,gy,gz) &
    !$acc& copyout(vortex,advection) if(gpu)
    call enforce_stairstep_compatibility_3d(rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,correction,imbalance)
    call solve_stairstep_poisson_3d(phi,rhs,solid,body,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
         nx,ny,nz,west,east,south,north,cart,tolerance,nsweep,merge(100,nsweep,tolerance>0),sweeps,residual,use_openacc=gpu, &
         method=trim(method))
    call MPI_Barrier(cart,ierr)
    t1=MPI_Wtime()
    call compute_stairstep_fpm_terms_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,solid,body,xc,yc,zc,dx,dy,dz, &
         xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vortex,advection,totals, &
         nx,ny,nz,west,east,south,north,cart,use_openacc=gpu)
    !$acc end data
    call MPI_Barrier(cart,ierr)
    t2=MPI_Wtime()
    call trace_pop()
    local_t=[t1-t0,t2-t1]
    call MPI_Allreduce(local_t,global_t,2,MPI_DOUBLE_PRECISION,MPI_MAX,cart,ierr)
    solve_seconds=global_t(1);terms_seconds=global_t(2)
    if(tolerance<=0.and.sweeps/=nsweep)call MPI_Abort(cart,4,ierr)
    if(tolerance>0.and.residual>tolerance)call MPI_Abort(cart,5,ierr)
  end subroutine run_case

  pure real(real64) function median3(values)
    real(real64),intent(in)::values(3)
    median3=sum(values)-minval(values)-maxval(values)
  end function median3
end program benchmark_stairstep_3d
