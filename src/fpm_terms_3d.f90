module fpm_terms_3d
  use,intrinsic :: iso_fortran_env,only:real64
  implicit none
  private
  public :: compute_fpm_terms_uniform_3d
  public :: compute_fpm_terms_nonuniform_3d
contains
  subroutine compute_fpm_terms_uniform_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,gx,gy,gz,vortex,advection,totals, &
       nx,ny,nz,dx,dy,dz)
    integer,intent(in)::nx,ny,nz,fluid(nx,ny,nz)
    real(real64),intent(in)::phi(0:nx+1,0:ny+1,0:nz+1),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
    real(real64),intent(in)::ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz),ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
    real(real64),intent(out)::gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz),vortex(nx,ny,nz),advection(nx,ny,nz),totals(2)
    real(real64),intent(in)::dx,dy,dz
    real(real64)::cx,cy,cz,mask,volume,vortex_total,advection_total
    integer::i,j,k
    vortex_total=0.0_real64;advection_total=0.0_real64;volume=dx*dy*dz
    !$acc data copyin(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid) copyout(gx,gy,gz,vortex,advection)
    !$acc parallel loop collapse(3) present(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,gx,gy,gz,vortex,advection) &
    !$acc& private(cx,cy,cz,mask) reduction(+:vortex_total,advection_total)
    do k=1,nz;do j=1,ny;do i=1,nx
      gx(i,j,k)=(phi(i+1,j,k)-phi(i-1,j,k))/(2.0_real64*dx)
      gy(i,j,k)=(phi(i,j+1,k)-phi(i,j-1,k))/(2.0_real64*dy)
      gz(i,j,k)=(phi(i,j,k+1)-phi(i,j,k-1))/(2.0_real64*dz)
      cx=w(i,j,k)*oy(i,j,k)-v(i,j,k)*oz(i,j,k)
      cy=u(i,j,k)*oz(i,j,k)-w(i,j,k)*ox(i,j,k)
      cz=v(i,j,k)*ox(i,j,k)-u(i,j,k)*oy(i,j,k)
      mask=real(fluid(i,j,k),real64)
      vortex(i,j,k)=-mask*(cx*gx(i,j,k)+cy*gy(i,j,k)+cz*gz(i,j,k))
      advection(i,j,k)=-mask*(ax(i,j,k)*gx(i,j,k)+ay(i,j,k)*gy(i,j,k)+az(i,j,k)*gz(i,j,k))
      vortex_total=vortex_total+volume*vortex(i,j,k)
      advection_total=advection_total+volume*advection(i,j,k)
    enddo;enddo;enddo
    !$acc end parallel loop
    !$acc end data
    totals=[vortex_total,advection_total]
  end subroutine compute_fpm_terms_uniform_3d

  subroutine compute_fpm_terms_nonuniform_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid, &
       xc,yc,zc,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi,bounds, &
       gx,gy,gz,vortex,advection,totals,nx,ny,nz)
    integer,intent(in)::nx,ny,nz,fluid(nx,ny,nz)
    real(real64),intent(in)::phi(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
    real(real64),intent(in)::ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz)
    real(real64),intent(in)::ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
    real(real64),intent(in)::xc(nx),yc(ny),zc(nz),dx(nx),dy(ny),dz(nz)
    real(real64),intent(in)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz)
    real(real64),intent(in)::zlo(nx,ny),zhi(nx,ny),bounds(6)
    real(real64),intent(out)::gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz)
    real(real64),intent(out)::vortex(nx,ny,nz),advection(nx,ny,nz),totals(2)
    real(real64)::cx,cy,cz,mask,inside,volume,vortex_total,advection_total
    real(real64)::west,east,south,north,bottom,top
    integer::i,j,k
    if(nx<2.or.ny<2.or.nz<2) error stop 'nonuniform 3-D FPM kernel requires at least two cells per direction'
    vortex_total=0.0_real64;advection_total=0.0_real64
    !$acc data copyin(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,xc,yc,zc,dx,dy,dz) &
    !$acc& copyin(xlo,xhi,ylo,yhi,zlo,zhi,bounds) copyout(gx,gy,gz,vortex,advection)
    !$acc parallel loop collapse(3) present(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,xc,yc,zc,dx,dy,dz) &
    !$acc& present(xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vortex,advection) &
    !$acc& private(cx,cy,cz,mask,inside,volume,west,east,south,north,bottom,top) &
    !$acc& reduction(+:vortex_total,advection_total)
    do k=1,nz;do j=1,ny;do i=1,nx
      if(i==1)then;west=xlo(j,k);else;west=(phi(i,j,k)-phi(i-1,j,k))/(xc(i)-xc(i-1));endif
      if(i==nx)then;east=xhi(j,k);else;east=(phi(i+1,j,k)-phi(i,j,k))/(xc(i+1)-xc(i));endif
      if(j==1)then;south=ylo(i,k);else;south=(phi(i,j,k)-phi(i,j-1,k))/(yc(j)-yc(j-1));endif
      if(j==ny)then;north=yhi(i,k);else;north=(phi(i,j+1,k)-phi(i,j,k))/(yc(j+1)-yc(j));endif
      if(k==1)then;bottom=zlo(i,j);else;bottom=(phi(i,j,k)-phi(i,j,k-1))/(zc(k)-zc(k-1));endif
      if(k==nz)then;top=zhi(i,j);else;top=(phi(i,j,k+1)-phi(i,j,k))/(zc(k+1)-zc(k));endif
      gx(i,j,k)=0.5_real64*(west+east)
      gy(i,j,k)=0.5_real64*(south+north)
      gz(i,j,k)=0.5_real64*(bottom+top)
      cx=w(i,j,k)*oy(i,j,k)-v(i,j,k)*oz(i,j,k)
      cy=u(i,j,k)*oz(i,j,k)-w(i,j,k)*ox(i,j,k)
      cz=v(i,j,k)*ox(i,j,k)-u(i,j,k)*oy(i,j,k)
      mask=real(fluid(i,j,k),real64)
      vortex(i,j,k)=-mask*(cx*gx(i,j,k)+cy*gy(i,j,k)+cz*gz(i,j,k))
      advection(i,j,k)=-mask*(ax(i,j,k)*gx(i,j,k)+ay(i,j,k)*gy(i,j,k)+az(i,j,k)*gz(i,j,k))
      inside=merge(1.0_real64,0.0_real64,xc(i)>bounds(1).and.xc(i)<bounds(2).and. &
           yc(j)>bounds(3).and.yc(j)<bounds(4).and.zc(k)>bounds(5).and.zc(k)<bounds(6))
      volume=inside*dx(i)*dy(j)*dz(k)
      vortex_total=vortex_total+volume*vortex(i,j,k)
      advection_total=advection_total+volume*advection(i,j,k)
    enddo;enddo;enddo
    !$acc end parallel loop
    !$acc end data
    totals=[vortex_total,advection_total]
  end subroutine compute_fpm_terms_nonuniform_3d
end module fpm_terms_3d
