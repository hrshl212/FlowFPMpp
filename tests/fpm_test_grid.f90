module fpm_test_grid
  ! Grid and decomposition helpers shared by the tests, benchmark, and demo.
  use,intrinsic :: iso_fortran_env,only:real64
  use mpi
  implicit none
  private
  public::stretched_axis,local_axis,cartesian_layout
contains
  ! Global axis on [lo,hi] with n cells whose widths vary smoothly by the
  ! relative amplitude (0 gives a uniform axis).
  subroutine stretched_axis(n,lo,hi,amplitude,face,center,width)
    integer,intent(in)::n
    real(real64),intent(in)::lo,hi,amplitude
    real(real64),intent(out)::face(0:n),center(n),width(n)
    real(real64)::weight(n)
    integer::m
    do m=1,n
      weight(m)=1+amplitude*sin(2*acos(-1.0_real64)*(m-0.5_real64)/n)
    enddo
    face(0)=lo
    do m=1,n
      width(m)=(hi-lo)*weight(m)/sum(weight)
      face(m)=face(m-1)+width(m)
      center(m)=0.5_real64*(face(m-1)+face(m))
    enddo
    face(n)=hi
  end subroutine stretched_axis

  ! Rank-local centres with one halo centre on each side (mirrored across a
  ! physical boundary) and local widths, for global cells offset+1:offset+n.
  subroutine local_axis(gn,gcenter,gwidth,offset,n,center,width)
    integer,intent(in)::gn,offset,n
    real(real64),intent(in)::gcenter(gn),gwidth(gn)
    real(real64),intent(out)::center(0:n+1),width(n)
    integer::m,g
    do m=0,n+1
      g=offset+m
      if(g<1)then;center(m)=gcenter(1)-gwidth(1)
      elseif(g>gn)then;center(m)=gcenter(gn)+gwidth(gn)
      else;center(m)=gcenter(g);endif
    enddo
    width=gwidth(offset+1:offset+n)
  end subroutine local_axis

  ! Creates a 2-D non-periodic Cartesian communicator over comm and returns the
  ! local block size, offsets, and neighbours. Aborts if the grid does not divide.
  subroutine cartesian_layout(comm,gnx,gny,cart,nx,ny,ioff,joff,west,east,south,north)
    integer,intent(in)::comm,gnx,gny
    integer,intent(out)::cart,nx,ny,ioff,joff,west,east,south,north
    integer::nproc,rank,dims(2),coords(2),ierr
    call MPI_Comm_size(comm,nproc,ierr)
    dims=0
    call MPI_Dims_create(nproc,2,dims,ierr)
    if(mod(gnx,dims(1))/=0.or.mod(gny,dims(2))/=0)call MPI_Abort(comm,2,ierr)
    call MPI_Cart_create(comm,2,dims,[.false.,.false.],.false.,cart,ierr)
    call MPI_Comm_rank(cart,rank,ierr)
    call MPI_Cart_coords(cart,rank,2,coords,ierr)
    call MPI_Cart_shift(cart,0,1,west,east,ierr)
    call MPI_Cart_shift(cart,1,1,south,north,ierr)
    nx=gnx/dims(1);ny=gny/dims(2);ioff=coords(1)*nx;joff=coords(2)*ny
  end subroutine cartesian_layout
end module fpm_test_grid
