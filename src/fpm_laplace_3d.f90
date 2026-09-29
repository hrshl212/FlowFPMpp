module fpm_laplace_3d
  use,intrinsic :: iso_fortran_env,only:real64
  implicit none
  private
  public :: enforce_neumann_compatibility_3d,solve_neumann_poisson_3d
  public :: enforce_neumann_compatibility_nonuniform_3d
  public :: solve_neumann_poisson_nonuniform_3d
contains
  subroutine enforce_neumann_compatibility_3d(rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz,correction,imbalance)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(in)::rhs(nx,ny,nz),dx,dy,dz
    real(real64),intent(inout)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(out)::correction,imbalance
    real(real64)::source,flux,area
    integer::i,j,k
    source=0.0_real64;flux=0.0_real64
    !$acc data copyin(rhs) copy(xlo,xhi,ylo,yhi,zlo,zhi)
    !$acc parallel loop collapse(3) present(rhs) reduction(+:source)
    do k=1,nz;do j=1,ny;do i=1,nx;source=source+rhs(i,j,k)*dx*dy*dz;enddo;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(xlo,xhi) reduction(+:flux)
    do k=1,nz;do j=1,ny;flux=flux+(-xlo(j,k)+xhi(j,k))*dy*dz;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi) reduction(+:flux)
    do k=1,nz;do i=1,nx;flux=flux+(-ylo(i,k)+yhi(i,k))*dx*dz;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi) reduction(+:flux)
    do j=1,ny;do i=1,nx;flux=flux+(-zlo(i,j)+zhi(i,j))*dx*dy;enddo;enddo
    !$acc end parallel loop
    area=2.0_real64*(real(ny*nz,real64)*dy*dz+real(nx*nz,real64)*dx*dz+real(nx*ny,real64)*dx*dy)
    correction=-(flux-source)/area
    !$acc parallel loop collapse(2) present(xlo,xhi)
    do k=1,nz;do j=1,ny;xlo(j,k)=xlo(j,k)-correction;xhi(j,k)=xhi(j,k)+correction;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi)
    do k=1,nz;do i=1,nx;ylo(i,k)=ylo(i,k)-correction;yhi(i,k)=yhi(i,k)+correction;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi)
    do j=1,ny;do i=1,nx;zlo(i,j)=zlo(i,j)-correction;zhi(i,j)=zhi(i,j)+correction;enddo;enddo
    !$acc end parallel loop
    flux=0.0_real64
    !$acc parallel loop collapse(2) present(xlo,xhi) reduction(+:flux)
    do k=1,nz;do j=1,ny;flux=flux+(-xlo(j,k)+xhi(j,k))*dy*dz;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi) reduction(+:flux)
    do k=1,nz;do i=1,nx;flux=flux+(-ylo(i,k)+yhi(i,k))*dx*dz;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi) reduction(+:flux)
    do j=1,ny;do i=1,nx;flux=flux+(-zlo(i,j)+zhi(i,j))*dx*dy;enddo;enddo
    !$acc end parallel loop
    imbalance=flux-source
    !$acc end data
  end subroutine enforce_neumann_compatibility_3d

  subroutine solve_neumann_poisson_3d(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz,tolerance, &
       max_iterations,check_interval,iterations,final_residual)
    integer,intent(in)::nx,ny,nz,max_iterations,check_interval
    real(real64),intent(inout)::phi(0:nx+1,0:ny+1,0:nz+1)
    real(real64),intent(in)::rhs(nx,ny,nz),xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(in)::dx,dy,dz,tolerance
    integer,intent(out)::iterations
    real(real64),intent(out)::final_residual
    real(real64),allocatable::next(:,:,:)
    real(real64)::idx2,idy2,idz2,invdiag,xm,xp,ym,yp,zm,zp,mean_value
    integer::i,j,k,step,nsteps
    allocate(next(0:nx+1,0:ny+1,0:nz+1));next=0.0_real64
    idx2=1.0_real64/dx**2;idy2=1.0_real64/dy**2;idz2=1.0_real64/dz**2
    invdiag=1.0_real64/(2.0_real64*(idx2+idy2+idz2))
    iterations=0;final_residual=huge(1.0_real64)
    !$acc data copy(phi) copyin(rhs,xlo,xhi,ylo,yhi,zlo,zhi) create(next)
    do while(iterations<max_iterations.and.final_residual>tolerance)
      nsteps=min(max(1,check_interval),max_iterations-iterations)
      do step=1,nsteps
        !$acc parallel loop collapse(3) present(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,next) private(xm,xp,ym,yp,zm,zp)
        do k=1,nz;do j=1,ny;do i=1,nx
          xm=phi(max(1,i-1),j,k);xp=phi(min(nx,i+1),j,k)
          ym=phi(i,max(1,j-1),k);yp=phi(i,min(ny,j+1),k)
          zm=phi(i,j,max(1,k-1));zp=phi(i,j,min(nz,k+1))
          if(i==1)xm=phi(i,j,k)-xlo(j,k)*dx
          if(i==nx)xp=phi(i,j,k)+xhi(j,k)*dx
          if(j==1)ym=phi(i,j,k)-ylo(i,k)*dy
          if(j==ny)yp=phi(i,j,k)+yhi(i,k)*dy
          if(k==1)zm=phi(i,j,k)-zlo(i,j)*dz
          if(k==nz)zp=phi(i,j,k)+zhi(i,j)*dz
          next(i,j,k)=((xm+xp)*idx2+(ym+yp)*idy2+(zm+zp)*idz2-rhs(i,j,k))*invdiag
        enddo;enddo;enddo
        !$acc end parallel loop
        !$acc parallel loop collapse(3) present(phi,next)
        do k=1,nz;do j=1,ny;do i=1,nx;phi(i,j,k)=next(i,j,k);enddo;enddo;enddo
        !$acc end parallel loop
        iterations=iterations+1
      enddo
      final_residual=residual_norm(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz)
    enddo
    mean_value=0.0_real64
    !$acc parallel loop collapse(3) present(phi) reduction(+:mean_value)
    do k=1,nz;do j=1,ny;do i=1,nx;mean_value=mean_value+phi(i,j,k);enddo;enddo;enddo
    !$acc end parallel loop
    mean_value=mean_value/real(nx*ny*nz,real64)
    !$acc parallel loop collapse(3) present(phi)
    do k=1,nz;do j=1,ny;do i=1,nx;phi(i,j,k)=phi(i,j,k)-mean_value;enddo;enddo;enddo
    !$acc end parallel loop
    call set_ghosts(phi,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz)
    final_residual=residual_norm(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz)
    !$acc end data
    deallocate(next)
  end subroutine solve_neumann_poisson_3d

  function residual_norm(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz) result(value)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny),dx,dy,dz
    real(real64)::value,xm,xp,ym,yp,zm,zp,res
    integer::i,j,k
    value=0.0_real64
    !$acc parallel loop collapse(3) present(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi) &
    !$acc& reduction(max:value) private(xm,xp,ym,yp,zm,zp,res)
    do k=1,nz;do j=1,ny;do i=1,nx
      xm=phi(max(1,i-1),j,k);xp=phi(min(nx,i+1),j,k)
      ym=phi(i,max(1,j-1),k);yp=phi(i,min(ny,j+1),k)
      zm=phi(i,j,max(1,k-1));zp=phi(i,j,min(nz,k+1))
      if(i==1)xm=phi(i,j,k)-xlo(j,k)*dx
      if(i==nx)xp=phi(i,j,k)+xhi(j,k)*dx
      if(j==1)ym=phi(i,j,k)-ylo(i,k)*dy
      if(j==ny)yp=phi(i,j,k)+yhi(i,k)*dy
      if(k==1)zm=phi(i,j,k)-zlo(i,j)*dz
      if(k==nz)zp=phi(i,j,k)+zhi(i,j)*dz
      res=(xm-2*phi(i,j,k)+xp)/dx**2+(ym-2*phi(i,j,k)+yp)/dy**2+(zm-2*phi(i,j,k)+zp)/dz**2-rhs(i,j,k)
      value=max(value,abs(res))
    enddo;enddo;enddo
    !$acc end parallel loop
  end function residual_norm

  subroutine set_ghosts(phi,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,dx,dy,dz)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(inout)::phi(0:nx+1,0:ny+1,0:nz+1)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny),dx,dy,dz
    integer::i,j,k
    !$acc parallel loop collapse(2) present(phi,xlo,xhi)
    do k=1,nz;do j=1,ny;phi(0,j,k)=phi(1,j,k)-xlo(j,k)*dx;phi(nx+1,j,k)=phi(nx,j,k)+xhi(j,k)*dx;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(phi,ylo,yhi)
    do k=1,nz;do i=1,nx;phi(i,0,k)=phi(i,1,k)-ylo(i,k)*dy;phi(i,ny+1,k)=phi(i,ny,k)+yhi(i,k)*dy;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(phi,zlo,zhi)
    do j=1,ny;do i=1,nx;phi(i,j,0)=phi(i,j,1)-zlo(i,j)*dz;phi(i,j,nz+1)=phi(i,j,nz)+zhi(i,j)*dz;enddo;enddo
    !$acc end parallel loop
  end subroutine set_ghosts

  subroutine enforce_neumann_compatibility_nonuniform_3d(rhs,xlo,xhi,ylo,yhi,zlo,zhi, &
       dx,dy,dz,nx,ny,nz,correction,imbalance)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(in)::rhs(nx,ny,nz),dx(nx),dy(ny),dz(nz)
    real(real64),intent(inout)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
    real(real64),intent(out)::correction,imbalance
    real(real64)::source,flux,area
    integer::i,j,k
    source=0.0_real64;flux=0.0_real64;area=0.0_real64
    !$acc data copyin(rhs,dx,dy,dz) copy(xlo,xhi,ylo,yhi,zlo,zhi)
    !$acc parallel loop collapse(3) present(rhs,dx,dy,dz) reduction(+:source)
    do k=1,nz;do j=1,ny;do i=1,nx
      source=source+rhs(i,j,k)*dx(i)*dy(j)*dz(k)
    enddo;enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(xlo,xhi,dy,dz) reduction(+:flux,area)
    do k=1,nz;do j=1,ny
      flux=flux+(-xlo(j,k)+xhi(j,k))*dy(j)*dz(k)
      area=area+2.0_real64*dy(j)*dz(k)
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi,dx,dz) reduction(+:flux,area)
    do k=1,nz;do i=1,nx
      flux=flux+(-ylo(i,k)+yhi(i,k))*dx(i)*dz(k)
      area=area+2.0_real64*dx(i)*dz(k)
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi,dx,dy) reduction(+:flux,area)
    do j=1,ny;do i=1,nx
      flux=flux+(-zlo(i,j)+zhi(i,j))*dx(i)*dy(j)
      area=area+2.0_real64*dx(i)*dy(j)
    enddo;enddo
    !$acc end parallel loop
    correction=-(flux-source)/area
    !$acc parallel loop collapse(2) present(xlo,xhi)
    do k=1,nz;do j=1,ny
      xlo(j,k)=xlo(j,k)-correction;xhi(j,k)=xhi(j,k)+correction
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi)
    do k=1,nz;do i=1,nx
      ylo(i,k)=ylo(i,k)-correction;yhi(i,k)=yhi(i,k)+correction
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi)
    do j=1,ny;do i=1,nx
      zlo(i,j)=zlo(i,j)-correction;zhi(i,j)=zhi(i,j)+correction
    enddo;enddo
    !$acc end parallel loop
    flux=0.0_real64
    !$acc parallel loop collapse(2) present(xlo,xhi,dy,dz) reduction(+:flux)
    do k=1,nz;do j=1,ny
      flux=flux+(-xlo(j,k)+xhi(j,k))*dy(j)*dz(k)
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(ylo,yhi,dx,dz) reduction(+:flux)
    do k=1,nz;do i=1,nx
      flux=flux+(-ylo(i,k)+yhi(i,k))*dx(i)*dz(k)
    enddo;enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(zlo,zhi,dx,dy) reduction(+:flux)
    do j=1,ny;do i=1,nx
      flux=flux+(-zlo(i,j)+zhi(i,j))*dx(i)*dy(j)
    enddo;enddo
    !$acc end parallel loop
    imbalance=flux-source
    !$acc end data
  end subroutine enforce_neumann_compatibility_nonuniform_3d

  subroutine solve_neumann_poisson_nonuniform_3d(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi, &
       xc,yc,zc,dx,dy,dz,nx,ny,nz,tolerance,max_iterations,check_interval,iterations,final_residual)
    integer,intent(in)::nx,ny,nz,max_iterations,check_interval
    real(real64),intent(inout)::phi(nx,ny,nz)
    real(real64),intent(in)::rhs(nx,ny,nz),xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz)
    real(real64),intent(in)::zlo(nx,ny),zhi(nx,ny),xc(nx),yc(ny),zc(nz),dx(nx),dy(ny),dz(nz),tolerance
    integer,intent(out)::iterations
    real(real64),intent(out)::final_residual
    real(real64),allocatable::next(:,:,:)
    real(real64)::aw,ae,as,an,ab,at,diag,numer,weighted_sum,volume_sum,mean_value
    integer::i,j,k,step,nsteps
    allocate(next(nx,ny,nz));next=0.0_real64
    iterations=0;final_residual=huge(1.0_real64)
    !$acc data copy(phi) copyin(rhs,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz) create(next)
    do while(iterations<max_iterations.and.final_residual>tolerance)
      nsteps=min(max(1,check_interval),max_iterations-iterations)
      do step=1,nsteps
        !$acc parallel loop collapse(3) present(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz,next) &
        !$acc& private(aw,ae,as,an,ab,at,diag,numer)
        do k=1,nz;do j=1,ny;do i=1,nx
          aw=0.0_real64;ae=0.0_real64;as=0.0_real64
          an=0.0_real64;ab=0.0_real64;at=0.0_real64
          if(i>1)aw=1.0_real64/(dx(i)*(xc(i)-xc(i-1)))
          if(i<nx)ae=1.0_real64/(dx(i)*(xc(i+1)-xc(i)))
          if(j>1)as=1.0_real64/(dy(j)*(yc(j)-yc(j-1)))
          if(j<ny)an=1.0_real64/(dy(j)*(yc(j+1)-yc(j)))
          if(k>1)ab=1.0_real64/(dz(k)*(zc(k)-zc(k-1)))
          if(k<nz)at=1.0_real64/(dz(k)*(zc(k+1)-zc(k)))
          diag=aw+ae+as+an+ab+at
          numer=-rhs(i,j,k)
          if(i>1)numer=numer+aw*phi(i-1,j,k)
          if(i<nx)numer=numer+ae*phi(i+1,j,k)
          if(j>1)numer=numer+as*phi(i,j-1,k)
          if(j<ny)numer=numer+an*phi(i,j+1,k)
          if(k>1)numer=numer+ab*phi(i,j,k-1)
          if(k<nz)numer=numer+at*phi(i,j,k+1)
          if(i==1)numer=numer-xlo(j,k)/dx(i)
          if(i==nx)numer=numer+xhi(j,k)/dx(i)
          if(j==1)numer=numer-ylo(i,k)/dy(j)
          if(j==ny)numer=numer+yhi(i,k)/dy(j)
          if(k==1)numer=numer-zlo(i,j)/dz(k)
          if(k==nz)numer=numer+zhi(i,j)/dz(k)
          ! Damping removes the checkerboard two-cycle of undamped Jacobi on
          ! the conservative Neumann operator's bipartite cell graph.
          next(i,j,k)=0.2_real64*phi(i,j,k)+0.8_real64*numer/diag
        enddo;enddo;enddo
        !$acc end parallel loop
        !$acc parallel loop collapse(3) present(phi,next)
        do k=1,nz;do j=1,ny;do i=1,nx
          phi(i,j,k)=next(i,j,k)
        enddo;enddo;enddo
        !$acc end parallel loop
        iterations=iterations+1
      enddo
      final_residual=residual_nonuniform(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi, &
           xc,yc,zc,dx,dy,dz,nx,ny,nz)
    enddo
    weighted_sum=0.0_real64;volume_sum=0.0_real64
    !$acc parallel loop collapse(3) present(phi,dx,dy,dz) reduction(+:weighted_sum,volume_sum)
    do k=1,nz;do j=1,ny;do i=1,nx
      weighted_sum=weighted_sum+phi(i,j,k)*dx(i)*dy(j)*dz(k)
      volume_sum=volume_sum+dx(i)*dy(j)*dz(k)
    enddo;enddo;enddo
    !$acc end parallel loop
    mean_value=weighted_sum/volume_sum
    !$acc parallel loop collapse(3) present(phi)
    do k=1,nz;do j=1,ny;do i=1,nx
      phi(i,j,k)=phi(i,j,k)-mean_value
    enddo;enddo;enddo
    !$acc end parallel loop
    final_residual=residual_nonuniform(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi, &
         xc,yc,zc,dx,dy,dz,nx,ny,nz)
    !$acc end data
    deallocate(next)
  end subroutine solve_neumann_poisson_nonuniform_3d

  function residual_nonuniform(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz,nx,ny,nz) result(value)
    integer,intent(in)::nx,ny,nz
    real(real64),intent(in)::phi(nx,ny,nz),rhs(nx,ny,nz),xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz)
    real(real64),intent(in)::zlo(nx,ny),zhi(nx,ny),xc(nx),yc(ny),zc(nz),dx(nx),dy(ny),dz(nz)
    real(real64)::value,gw,ge,gs,gn,gb,gt,res
    integer::i,j,k
    value=0.0_real64
    !$acc parallel loop collapse(3) present(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz) &
    !$acc& private(gw,ge,gs,gn,gb,gt,res) reduction(max:value)
    do k=1,nz;do j=1,ny;do i=1,nx
      if(i==1)then;gw=xlo(j,k);else;gw=(phi(i,j,k)-phi(i-1,j,k))/(xc(i)-xc(i-1));endif
      if(i==nx)then;ge=xhi(j,k);else;ge=(phi(i+1,j,k)-phi(i,j,k))/(xc(i+1)-xc(i));endif
      if(j==1)then;gs=ylo(i,k);else;gs=(phi(i,j,k)-phi(i,j-1,k))/(yc(j)-yc(j-1));endif
      if(j==ny)then;gn=yhi(i,k);else;gn=(phi(i,j+1,k)-phi(i,j,k))/(yc(j+1)-yc(j));endif
      if(k==1)then;gb=zlo(i,j);else;gb=(phi(i,j,k)-phi(i,j,k-1))/(zc(k)-zc(k-1));endif
      if(k==nz)then;gt=zhi(i,j);else;gt=(phi(i,j,k+1)-phi(i,j,k))/(zc(k+1)-zc(k));endif
      res=(ge-gw)/dx(i)+(gn-gs)/dy(j)+(gt-gb)/dz(k)-rhs(i,j,k)
      value=max(value,abs(res))
    enddo;enddo;enddo
    !$acc end parallel loop
  end function residual_nonuniform
end module fpm_laplace_3d
