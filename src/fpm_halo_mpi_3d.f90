module fpm_halo_mpi_3d
  use,intrinsic :: iso_fortran_env,only:real64
  use mpi
  implicit none
  private
  type,public::halo_workspace_3d
    integer::nx=0,ny=0,nz=0,depth=0
    logical::device_x=.false.,temporary_x=.false.,unbatched_x=.false.
    real(real64),allocatable::send_w(:),send_e(:),recv_w(:),recv_e(:)
    real(real64),allocatable::send_s(:),send_n(:),recv_s(:),recv_n(:)
  end type halo_workspace_3d
  type,public::integer_halo_workspace_3d
    integer::nx=0,ny=0,nz=0,depth=0
    integer,allocatable::send_w(:),send_e(:),recv_w(:),recv_e(:)
    integer,allocatable::send_s(:),send_n(:),recv_s(:),recv_n(:)
  end type integer_halo_workspace_3d
  public::initialize_halo_workspace_3d,finalize_halo_workspace_3d,exchange_xy_halo_3d
  public::initialize_integer_halo_workspace_3d,finalize_integer_halo_workspace_3d,exchange_integer_xy_halo_3d
  public::exchange_integer_xy_halo_padded_3d
contains
  subroutine initialize_integer_halo_workspace_3d(work,nx,ny,nz,depth)
    type(integer_halo_workspace_3d),intent(inout)::work
    integer,intent(in)::nx,ny,nz,depth
    call finalize_integer_halo_workspace_3d(work)
    if(depth<1.or.depth>min(nx,ny))error stop 'invalid integer 3-D halo depth'
    work%nx=nx;work%ny=ny;work%nz=nz;work%depth=depth
    allocate(work%send_w(depth*ny*nz),work%send_e(depth*ny*nz),work%recv_w(depth*ny*nz),work%recv_e(depth*ny*nz))
    allocate(work%send_s((nx+2*depth)*depth*nz),work%send_n((nx+2*depth)*depth*nz), &
         work%recv_s((nx+2*depth)*depth*nz),work%recv_n((nx+2*depth)*depth*nz))
  end subroutine initialize_integer_halo_workspace_3d

  subroutine finalize_integer_halo_workspace_3d(work)
    type(integer_halo_workspace_3d),intent(inout)::work
    if(allocated(work%send_w))deallocate(work%send_w,work%send_e,work%recv_w,work%recv_e, &
         work%send_s,work%send_n,work%recv_s,work%recv_n)
    work%nx=0;work%ny=0;work%nz=0;work%depth=0
  end subroutine finalize_integer_halo_workspace_3d

  subroutine exchange_integer_xy_halo_3d(field,nx,ny,nz,depth,west,east,south,north,comm,work)
    integer,intent(in)::nx,ny,nz,depth,west,east,south,north,comm
    integer,intent(inout)::field(1-depth:nx+depth,1-depth:ny+depth,nz)
    type(integer_halo_workspace_3d),intent(inout)::work
    call exchange_integer_xy_halo_padded_3d(field,nx,ny,nz,depth,depth,west,east,south,north,comm,work)
  end subroutine exchange_integer_xy_halo_3d

  ! Exchanges depth layers of a field padded by pad>=depth layers, without
  ! requiring callers to pass a (copied) non-contiguous depth-sized section.
  subroutine exchange_integer_xy_halo_padded_3d(field,nx,ny,nz,pad,depth,west,east,south,north,comm,work)
    integer,intent(in)::nx,ny,nz,pad,depth,west,east,south,north,comm
    integer,intent(inout)::field(1-pad:nx+pad,1-pad:ny+pad,nz)
    type(integer_halo_workspace_3d),intent(inout)::work
    integer::i,j,k,q,ierr,status(MPI_STATUS_SIZE),xcount,ycount
    if(any([work%nx,work%ny,work%nz,work%depth]/=[nx,ny,nz,depth]))error stop 'integer halo workspace mismatch'
    if(pad<depth)error stop 'integer halo padding smaller than exchange depth'
    if(all([west,east,south,north]==MPI_PROC_NULL))return
    q=0
    do k=1,nz;do j=1,ny;do i=1,depth
      q=q+1;work%send_w(q)=field(i,j,k);work%send_e(q)=field(nx-depth+i,j,k)
    enddo;enddo;enddo;xcount=q
    call MPI_Sendrecv(work%send_w,xcount,MPI_INTEGER,west,411,work%recv_e,xcount,MPI_INTEGER,east,411,comm,status,ierr)
    call MPI_Sendrecv(work%send_e,xcount,MPI_INTEGER,east,412,work%recv_w,xcount,MPI_INTEGER,west,412,comm,status,ierr)
    q=0
    do k=1,nz;do j=1,ny;do i=1,depth
      q=q+1
      if(west/=MPI_PROC_NULL)field(i-depth,j,k)=work%recv_w(q)
      if(east/=MPI_PROC_NULL)field(nx+i,j,k)=work%recv_e(q)
    enddo;enddo;enddo
    q=0
    do k=1,nz;do j=1,depth;do i=1-depth,nx+depth
      q=q+1;work%send_s(q)=field(i,j,k);work%send_n(q)=field(i,ny-depth+j,k)
    enddo;enddo;enddo;ycount=q
    call MPI_Sendrecv(work%send_s,ycount,MPI_INTEGER,south,413,work%recv_n,ycount,MPI_INTEGER,north,413,comm,status,ierr)
    call MPI_Sendrecv(work%send_n,ycount,MPI_INTEGER,north,414,work%recv_s,ycount,MPI_INTEGER,south,414,comm,status,ierr)
    q=0
    do k=1,nz;do j=1,depth;do i=1-depth,nx+depth
      q=q+1
      if(south/=MPI_PROC_NULL)field(i,j-depth,k)=work%recv_s(q)
      if(north/=MPI_PROC_NULL)field(i,ny+j,k)=work%recv_n(q)
    enddo;enddo;enddo
  end subroutine exchange_integer_xy_halo_padded_3d

  subroutine initialize_halo_workspace_3d(work,nx,ny,nz,depth)
    type(halo_workspace_3d),intent(inout)::work
    integer,intent(in)::nx,ny,nz,depth
    character(len=8)::setting
    integer::length,status
    call finalize_halo_workspace_3d(work)
    if(depth<1.or.depth>min(nx,ny))error stop 'invalid 3-D halo depth'
    work%nx=nx;work%ny=ny;work%nz=nz;work%depth=depth
    allocate(work%send_w(depth*ny*nz),work%send_e(depth*ny*nz),work%recv_w(depth*ny*nz),work%recv_e(depth*ny*nz))
    allocate(work%send_s((nx+2*depth)*depth*nz),work%send_n((nx+2*depth)*depth*nz), &
         work%recv_s((nx+2*depth)*depth*nz),work%recv_n((nx+2*depth)*depth*nz))
    ! Benchmark-only switch restoring per-exchange temporary x-face device buffers.
    call get_environment_variable('FPM_HALO_TEMPORARY_BUFFERS',setting,length,status)
    work%temporary_x=status==0.and.length>0.and.setting(1:1)=='1'
    ! Benchmark-only switch restoring the pre-batching per-directive x-face updates.
    call get_environment_variable('FPM_HALO_UNBATCHED_X',setting,length,status)
    work%unbatched_x=status==0.and.length>0.and.setting(1:1)=='1'
  end subroutine initialize_halo_workspace_3d

  subroutine finalize_halo_workspace_3d(work)
    type(halo_workspace_3d),intent(inout)::work
    if(work%device_x)call delete_device_x_faces(work%send_w,work%send_e,work%recv_w,work%recv_e,size(work%send_w))
    work%device_x=.false.;work%temporary_x=.false.;work%unbatched_x=.false.
    if(allocated(work%send_w))deallocate(work%send_w,work%send_e,work%recv_w,work%recv_e, &
         work%send_s,work%send_n,work%recv_s,work%recv_n)
    work%nx=0;work%ny=0;work%nz=0;work%depth=0
  end subroutine finalize_halo_workspace_3d

  subroutine create_device_x_faces(send_w,send_e,recv_w,recv_e,n)
    integer,intent(in)::n
    real(real64),intent(inout)::send_w(n),send_e(n),recv_w(n),recv_e(n)
    !$acc enter data create(send_w,send_e,recv_w,recv_e)
  end subroutine create_device_x_faces

  subroutine delete_device_x_faces(send_w,send_e,recv_w,recv_e,n)
    integer,intent(in)::n
    real(real64),intent(inout)::send_w(n),send_e(n),recv_w(n),recv_e(n)
    !$acc exit data delete(send_w,send_e,recv_w,recv_e)
  end subroutine delete_device_x_faces

  subroutine exchange_xy_halo_3d(field,nx,ny,nz,depth,west,east,south,north,comm,work,host_only,batched_device)
    integer,intent(in)::nx,ny,nz,depth,west,east,south,north,comm
    real(real64),intent(inout)::field(1-depth:nx+depth,1-depth:ny+depth,nz)
    type(halo_workspace_3d),intent(inout)::work
    integer::i,j,k,q,ierr,status(MPI_STATUS_SIZE),xcount,ycount
    real(real64),allocatable::packed_send_w(:),packed_send_e(:),packed_recv_w(:),packed_recv_e(:)
    logical,intent(in),optional::host_only,batched_device
    logical::skip_device,use_batch
    if(any([work%nx,work%ny,work%nz,work%depth]/=[nx,ny,nz,depth]))error stop '3-D halo workspace mismatch'
    skip_device=.false.;if(present(host_only))skip_device=host_only
    use_batch=.false.;if(present(batched_device))use_batch=batched_device.and..not.skip_device.and..not.work%unbatched_x
    if(all([west,east,south,north]==MPI_PROC_NULL))return
    xcount=depth*ny*nz
    if(use_batch.and.work%temporary_x)then
      allocate(packed_send_w(xcount),packed_send_e(xcount),packed_recv_w(xcount),packed_recv_e(xcount))
      !$acc data create(packed_send_w,packed_send_e,packed_recv_w,packed_recv_e)
      call exchange_x_faces_device(field,nx,ny,nz,depth,west,east,comm,packed_send_w,packed_send_e, &
           packed_recv_w,packed_recv_e,xcount)
      !$acc end data
      deallocate(packed_send_w,packed_send_e,packed_recv_w,packed_recv_e)
    elseif(use_batch)then
      if(.not.work%device_x)then
        call create_device_x_faces(work%send_w,work%send_e,work%recv_w,work%recv_e,xcount);work%device_x=.true.
      endif
      call exchange_x_faces_device(field,nx,ny,nz,depth,west,east,comm,work%send_w,work%send_e, &
           work%recv_w,work%recv_e,xcount)
    else
      !$acc update self(field(1:depth,1:ny,1:nz)) if_present if(west/=MPI_PROC_NULL.and..not.skip_device)
      !$acc update self(field(nx-depth+1:nx,1:ny,1:nz)) if_present if(east/=MPI_PROC_NULL.and..not.skip_device)
      q=0
      do k=1,nz;do j=1,ny;do i=1,depth
        q=q+1;work%send_w(q)=field(i,j,k);work%send_e(q)=field(nx-depth+i,j,k)
      enddo;enddo;enddo
      call MPI_Sendrecv(work%send_w,xcount,MPI_DOUBLE_PRECISION,west,401,work%recv_e,xcount,MPI_DOUBLE_PRECISION, &
           east,401,comm,status,ierr)
      call MPI_Sendrecv(work%send_e,xcount,MPI_DOUBLE_PRECISION,east,402,work%recv_w,xcount,MPI_DOUBLE_PRECISION, &
           west,402,comm,status,ierr)
      q=0
      do k=1,nz;do j=1,ny;do i=1,depth
        q=q+1
        if(west/=MPI_PROC_NULL)field(i-depth,j,k)=work%recv_w(q)
        if(east/=MPI_PROC_NULL)field(nx+i,j,k)=work%recv_e(q)
      enddo;enddo;enddo
      !$acc update device(field(1-depth:0,1:ny,1:nz)) if_present if(west/=MPI_PROC_NULL.and..not.skip_device)
      !$acc update device(field(nx+1:nx+depth,1:ny,1:nz)) if_present if(east/=MPI_PROC_NULL.and..not.skip_device)
    endif
    !$acc update self(field(1-depth:nx+depth,1:depth,1:nz)) if_present if(south/=MPI_PROC_NULL.and..not.skip_device)
    !$acc update self(field(1-depth:nx+depth,ny-depth+1:ny,1:nz)) if_present if(north/=MPI_PROC_NULL.and..not.skip_device)
    q=0
    do k=1,nz;do j=1,depth;do i=1-depth,nx+depth
      q=q+1;work%send_s(q)=field(i,j,k);work%send_n(q)=field(i,ny-depth+j,k)
    enddo;enddo;enddo
    ycount=q
    call MPI_Sendrecv(work%send_s,ycount,MPI_DOUBLE_PRECISION,south,403,work%recv_n,ycount,MPI_DOUBLE_PRECISION, &
         north,403,comm,status,ierr)
    call MPI_Sendrecv(work%send_n,ycount,MPI_DOUBLE_PRECISION,north,404,work%recv_s,ycount,MPI_DOUBLE_PRECISION, &
         south,404,comm,status,ierr)
    q=0
    do k=1,nz;do j=1,depth;do i=1-depth,nx+depth
      q=q+1
      if(south/=MPI_PROC_NULL)field(i,j-depth,k)=work%recv_s(q)
      if(north/=MPI_PROC_NULL)field(i,ny+j,k)=work%recv_n(q)
    enddo;enddo;enddo
    !$acc update device(field(1-depth:nx+depth,1-depth:0,1:nz)) if_present if(south/=MPI_PROC_NULL.and..not.skip_device)
    !$acc update device(field(1-depth:nx+depth,ny+1:ny+depth,1:nz)) if_present if(north/=MPI_PROC_NULL.and..not.skip_device)
  end subroutine exchange_xy_halo_3d

  subroutine exchange_x_faces_device(field,nx,ny,nz,depth,west,east,comm,packed_send_w,packed_send_e, &
       packed_recv_w,packed_recv_e,xcount)
    integer,intent(in)::nx,ny,nz,depth,west,east,comm,xcount
    real(real64),intent(inout)::field(1-depth:nx+depth,1-depth:ny+depth,nz)
    real(real64),intent(inout)::packed_send_w(xcount),packed_send_e(xcount),packed_recv_w(xcount),packed_recv_e(xcount)
    integer::i,j,k,q,ierr,status(MPI_STATUS_SIZE)
    !$acc parallel loop collapse(3) present(field,packed_send_w,packed_send_e) private(q)
    do k=1,nz;do j=1,ny;do i=1,depth
      q=((k-1)*ny+j-1)*depth+i
      if(west/=MPI_PROC_NULL)packed_send_w(q)=field(i,j,k)
      if(east/=MPI_PROC_NULL)packed_send_e(q)=field(nx-depth+i,j,k)
    enddo;enddo;enddo
    !$acc end parallel loop
    !$acc update self(packed_send_w) if(west/=MPI_PROC_NULL)
    !$acc update self(packed_send_e) if(east/=MPI_PROC_NULL)
    call MPI_Sendrecv(packed_send_w,xcount,MPI_DOUBLE_PRECISION,west,401,packed_recv_e,xcount, &
         MPI_DOUBLE_PRECISION,east,401,comm,status,ierr)
    call MPI_Sendrecv(packed_send_e,xcount,MPI_DOUBLE_PRECISION,east,402,packed_recv_w,xcount, &
         MPI_DOUBLE_PRECISION,west,402,comm,status,ierr)
    q=0
    do k=1,nz;do j=1,ny;do i=1,depth
      q=q+1
      if(west/=MPI_PROC_NULL)field(i-depth,j,k)=packed_recv_w(q)
      if(east/=MPI_PROC_NULL)field(nx+i,j,k)=packed_recv_e(q)
    enddo;enddo;enddo
    !$acc update device(packed_recv_w) if(west/=MPI_PROC_NULL)
    !$acc update device(packed_recv_e) if(east/=MPI_PROC_NULL)
    !$acc parallel loop collapse(3) present(field,packed_recv_w,packed_recv_e) private(q)
    do k=1,nz;do j=1,ny;do i=1,depth
      q=((k-1)*ny+j-1)*depth+i
      if(west/=MPI_PROC_NULL)field(i-depth,j,k)=packed_recv_w(q)
      if(east/=MPI_PROC_NULL)field(nx+i,j,k)=packed_recv_e(q)
    enddo;enddo;enddo
    !$acc end parallel loop
  end subroutine exchange_x_faces_device
end module fpm_halo_mpi_3d
