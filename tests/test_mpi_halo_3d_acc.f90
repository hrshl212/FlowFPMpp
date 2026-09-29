program test_mpi_halo_3d_acc
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use mpi
  use fpm_halo_mpi_3d
  implicit none
  integer,parameter::nx=5,ny=6,nz=3,depth=3
  real(real64),parameter::sentinel=-9999.0_real64
  real(real64)::field(1-depth:nx+depth,1-depth:ny+depth,nz),expected,local_error,error
  integer::rank,nproc,cart,coords(2),dims(2),west,east,south,north,i,j,k,gx,gy,ierr
  logical::periods(2),reorder
  type(halo_workspace_3d)::work

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=4)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  dims=[2,2];periods=.false.;reorder=.false.;call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,cart,ierr)
  call MPI_Cart_coords(cart,rank,2,coords,ierr);call MPI_Cart_shift(cart,0,1,west,east,ierr)
  call MPI_Cart_shift(cart,1,1,south,north,ierr)
  field=sentinel
  do k=1,nz;do j=1,ny;do i=1,nx
    gx=coords(1)*nx+i;gy=coords(2)*ny+j;field(i,j,k)=cell_value(gx,gy,k)
  enddo;enddo;enddo
  call initialize_halo_workspace_3d(work,nx,ny,nz,depth)

  !$acc data copy(field)
  call exchange_xy_halo_3d(field,nx,ny,nz,depth,west,east,south,north,cart,work,batched_device=.true.)
  local_error=0
  !$acc parallel loop collapse(3) present(field) reduction(max:local_error) private(gx,gy,expected)
  do k=1,nz;do j=1-depth,ny+depth;do i=1-depth,nx+depth
    gx=coords(1)*nx+i;gy=coords(2)*ny+j
    if(gx>=1.and.gx<=dims(1)*nx.and.gy>=1.and.gy<=dims(2)*ny)then
      expected=real(gx+100*gy+10000*k,real64)
    else
      expected=sentinel
    endif
    local_error=max(local_error,abs(field(i,j,k)-expected))
  enddo;enddo;enddo
  !$acc end parallel loop
  !$acc end data

  call MPI_Allreduce(local_error,error,1,MPI_DOUBLE_PRECISION,MPI_MAX,cart,ierr)
  if(rank==0)then
    write(output_unit,'(a,es14.6)')'OpenACC 3-D Cartesian device-halo error: ',error
    if(error==0)then;write(output_unit,'(a)')'PASS: post-MPI halo values are current on device'
    else;write(error_unit,'(a)')'FAIL: stale post-MPI device halo values';endif
  endif
  call finalize_halo_workspace_3d(work);call MPI_Comm_free(cart,ierr);call MPI_Finalize(ierr)
  if(error/=0)error stop 1
contains
  pure real(real64) function cell_value(ii,jj,kk)
    integer,intent(in)::ii,jj,kk
    cell_value=real(ii+100*jj+10000*kk,real64)
  end function cell_value
end program test_mpi_halo_3d_acc
