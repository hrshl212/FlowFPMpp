module fpm_partition_2d
  ! Force-partitioning terms (Menon & Mittal; Raut, Seo & Mittal 2024, eq. 2.7)
  ! for a 2-D stair-step body with rigid in-plane motion, per unit span and
  ! unit density, on one rank (host code):
  !   F_KINE  = - sum_B  (n . a_B) phi dS          body acceleration a_B
  !   F_VIF   =   sum_V  phi div[(u.grad)u] dV     vortex-induced force
  !   F_VIS   =   nu sum_V lap(u) . grad(phi) dV
  !   F_BOUND = - sum_S  (n . Du/Dt) phi dS        outer boundary
  ! n points out of the fluid. The paper writes F_VIF = -2 int Q phi dV; for
  ! incompressible flow div[(u.grad)u] = -2Q, and this divergence (face-flux)
  ! form is used because the pointwise Q of the discrete velocity is not the
  ! discrete divergence next to the wall, where phi is largest.
  !
  ! Discretization, on the same faces as the stair-step solver:
  ! - body faces (fluid/solid) carry the rigid-body velocity
  !   U_B = (k1 - k3 y, k2 + k3 x), acceleration
  !   a_B = (k4 - k6 y - k3^2 x, k5 + k6 x - k3^2 y), and convective
  !   acceleration (U_B.grad)U_B = k3 z x U_B; phi is extrapolated to the face
  !   with the prescribed face gradient;
  ! - outer faces use zero normal gradient for u and phi;
  ! - Du/Dt on the outer boundary is (u(t+dt) - u(t-dt))/(2 dt) + (u.grad)u.
  use,intrinsic :: iso_fortran_env,only:real64
  implicit none
  private
  public::fpm_partition_terms_2d
contains
  ! phi: influence field (fluid-volume gauge) for the direction whose body-face
  ! gradient is g=(gx,gy). ring(:,c,t): outer-boundary cells (west column,
  ! east column, south row, north row) of component c=1 (u), 2 (v) at t=1
  ! (previous) and t=2 (next snapshot). terms = [KINE, VIF, VIS, BOUND].
  ! Outputs the F_VIF density phi*div[(u.grad)u] and the vorticity (zero in solid).
  subroutine fpm_partition_terms_2d(phi,g,u,v,solid,xf,yf,kin,ring,dt,nu,nx,ny,terms,vif_density,omega)
    integer,intent(in)::nx,ny
    real(real64),intent(in)::phi(nx,ny),g(2),u(nx,ny),v(nx,ny),xf(0:nx),yf(0:ny),kin(6)
    real(real64),intent(in)::ring(2*nx+2*ny,2,2),dt,nu
    logical,intent(in)::solid(nx,ny)
    real(real64),intent(out)::terms(4),vif_density(nx,ny),omega(nx,ny)
    real(real64),allocatable::xc(:),yc(:),dx(:),dy(:),nxf(:,:),nyf(:,:),gpx(:,:),gpy(:,:)
    real(real64)::fu(4),fv(4),fp(4),div,lapu,lapv,hw,d,fx,fy,phif,ax,ay,ubx,uby,len,sgn,dudt,dvdt,nn
    integer::i,j,f,ii,jj,axis,m
    allocate(xc(nx),yc(ny),dx(nx),dy(ny),nxf(nx,ny),nyf(nx,ny),gpx(nx,ny),gpy(nx,ny))
    dx=xf(1:nx)-xf(0:nx-1);dy=yf(1:ny)-yf(0:ny-1)
    xc=0.5_real64*(xf(0:nx-1)+xf(1:nx));yc=0.5_real64*(yf(0:ny-1)+yf(1:ny))
    terms=0;vif_density=0;omega=0;nxf=0;nyf=0;gpx=0;gpy=0
    ! Pass 1: velocity gradients, (u.grad)u, vorticity, and grad(phi) per fluid cell.
    do j=1,ny;do i=1,nx
      if(solid(i,j))cycle
      do f=1,4
        call face(i,j,f,ii,jj,axis,sgn,hw,len,d,fx,fy)
        if(ii<1.or.ii>nx.or.jj<1.or.jj>ny)then          ! outer face
          fu(f)=u(i,j);fv(f)=v(i,j);fp(f)=0
        elseif(solid(ii,jj))then                       ! body face
          call body_velocity(fx,fy,fu(f),fv(f))
          fp(f)=g(axis)
        else                                           ! interior face, linear interpolation
          fu(f)=u(i,j)+(u(ii,jj)-u(i,j))*hw/d
          fv(f)=v(i,j)+(v(ii,jj)-v(i,j))*hw/d
          fp(f)=sgn*(phi(ii,jj)-phi(i,j))/d
        endif
      enddo
      ! faces: 1 east, 2 west, 3 north, 4 south
      nxf(i,j)=u(i,j)*(fu(1)-fu(2))/dx(i)+v(i,j)*(fu(3)-fu(4))/dy(j)
      nyf(i,j)=u(i,j)*(fv(1)-fv(2))/dx(i)+v(i,j)*(fv(3)-fv(4))/dy(j)
      omega(i,j)=(fv(1)-fv(2))/dx(i)-(fu(3)-fu(4))/dy(j)
      gpx(i,j)=0.5_real64*(fp(1)+fp(2));gpy(i,j)=0.5_real64*(fp(3)+fp(4))
    enddo;enddo
    ! Pass 2: the four terms.
    m=0
    do j=1,ny;do i=1,nx
      if(solid(i,j))cycle
      div=0;lapu=0;lapv=0
      do f=1,4
        call face(i,j,f,ii,jj,axis,sgn,hw,len,d,fx,fy)
        if(ii<1.or.ii>nx.or.jj<1.or.jj>ny)then
          ! Outer face: (u.grad)u from the cell, zero-gradient viscous flux, F_BOUND.
          nn=merge(nxf(i,j),nyf(i,j),axis==1)
          div=div+sgn*nn*len
          call ring_index(i,j,f,m)
          dudt=(ring(m,1,2)-ring(m,1,1))/(2*dt);dvdt=(ring(m,2,2)-ring(m,2,1))/(2*dt)
          terms(4)=terms(4)-sgn*merge(dudt+nxf(i,j),dvdt+nyf(i,j),axis==1)*phi(i,j)*len
        elseif(solid(ii,jj))then
          ! Body face: F_KINE, rigid convective acceleration, no-slip viscous flux.
          call body_velocity(fx,fy,ubx,uby)
          call body_acceleration(fx,fy,ax,ay)
          phif=phi(i,j)+sgn*hw*g(axis)
          terms(1)=terms(1)-sgn*merge(ax,ay,axis==1)*phif*len
          div=div+sgn*merge(-kin(3)*uby,kin(3)*ubx,axis==1)*len
          lapu=lapu+(ubx-u(i,j))/hw*len;lapv=lapv+(uby-v(i,j))/hw*len
        else
          div=div+sgn*merge(nxf(i,j)+(nxf(ii,jj)-nxf(i,j))*hw/d,nyf(i,j)+(nyf(ii,jj)-nyf(i,j))*hw/d,axis==1)*len
          lapu=lapu+(u(ii,jj)-u(i,j))/d*len;lapv=lapv+(v(ii,jj)-v(i,j))/d*len
        endif
      enddo
      vif_density(i,j)=phi(i,j)*div/(dx(i)*dy(j))
      terms(2)=terms(2)+phi(i,j)*div
      terms(3)=terms(3)+nu*(lapu*gpx(i,j)+lapv*gpy(i,j))
    enddo;enddo
    deallocate(xc,yc,dx,dy,nxf,nyf,gpx,gpy)
  contains
    ! Face f of cell (i,j): neighbour, axis, outward sign, half width, length,
    ! centre distance, and face-centre coordinates.
    subroutine face(i,j,f,ii,jj,axis,sgn,hw,len,d,fx,fy)
      integer,intent(in)::i,j,f
      integer,intent(out)::ii,jj,axis
      real(real64),intent(out)::sgn,hw,len,d,fx,fy
      select case(f)
      case(1);ii=i+1;jj=j;axis=1;sgn=1
      case(2);ii=i-1;jj=j;axis=1;sgn=-1
      case(3);ii=i;jj=j+1;axis=2;sgn=1
      case default;ii=i;jj=j-1;axis=2;sgn=-1
      end select
      if(axis==1)then
        hw=0.5_real64*dx(i);len=dy(j);fx=xc(i)+sgn*hw;fy=yc(j)
        d=hw;if(ii>=1.and.ii<=nx)d=abs(xc(ii)-xc(i))
      else
        hw=0.5_real64*dy(j);len=dx(i);fx=xc(i);fy=yc(j)+sgn*hw
        d=hw;if(jj>=1.and.jj<=ny)d=abs(yc(jj)-yc(j))
      endif
    end subroutine face

    subroutine ring_index(i,j,f,m)
      integer,intent(in)::i,j,f
      integer,intent(out)::m
      select case(f)
      case(2);m=j                  ! west column
      case(1);m=ny+j               ! east column
      case(4);m=2*ny+i             ! south row
      case default;m=2*ny+nx+i     ! north row
      end select
    end subroutine ring_index

    subroutine body_velocity(x,y,bu,bv)
      real(real64),intent(in)::x,y
      real(real64),intent(out)::bu,bv
      bu=kin(1)-kin(3)*y;bv=kin(2)+kin(3)*x
    end subroutine body_velocity

    subroutine body_acceleration(x,y,bx,by)
      real(real64),intent(in)::x,y
      real(real64),intent(out)::bx,by
      bx=kin(4)-kin(6)*y-kin(3)**2*x;by=kin(5)+kin(6)*x-kin(3)**2*y
    end subroutine body_acceleration
  end subroutine fpm_partition_terms_2d
end module fpm_partition_2d
