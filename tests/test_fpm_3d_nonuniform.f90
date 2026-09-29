program test_fpm_3d_nonuniform
  use,intrinsic :: iso_fortran_env,only:real64,output_unit,error_unit
  use fpm_laplace_3d
  use fpm_terms_3d
  implicit none
  integer,parameter::nx=13,ny=11,nz=9
  real(real64),parameter::tol=1.0e-10_real64
  real(real64)::xf(0:nx),yf(0:ny),zf(0:nz),xc(nx),yc(ny),zc(nz),dx(nx),dy(ny),dz(nz)
  real(real64)::phi(nx,ny,nz),exact(nx,ny,nz),rhs(nx,ny,nz)
  real(real64)::u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),ox(nx,ny,nz),oy(nx,ny,nz),oz(nx,ny,nz)
  real(real64)::ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
  real(real64)::gx(nx,ny,nz),gy(nx,ny,nz),gz(nx,ny,nz),vf(nx,ny,nz),af(nx,ny,nz)
  real(real64)::xlo(ny,nz),xhi(ny,nz),ylo(nx,nz),yhi(nx,nz),zlo(nx,ny),zhi(nx,ny)
  real(real64)::txlo(ny,nz),txhi(ny,nz),tylo(nx,nz),tyhi(nx,nz),tzlo(nx,ny),tzhi(nx,ny)
  integer::fluid(nx,ny,nz),i,j,k,iterations
  real(real64)::gw,ge,gs,gn,gb,gt,volume,total_volume,mean_exact
  real(real64)::correction,imbalance,residual,solution_error,gradient_error,field_error,integral_error
  real(real64)::egx,egy,egz,cx,cy,cz,ev,ea,expected_totals(2),totals(2),bounds(6),inside
  logical::passed

  call make_axis(xf,xc,dx,nx,0.23_real64)
  call make_axis(yf,yc,dy,ny,-0.17_real64)
  call make_axis(zf,zc,dz,nz,0.19_real64)
  do k=1,nz;do j=1,ny;do i=1,nx
    exact(i,j,k)=field(xc(i),yc(j),zc(k))
  enddo;enddo;enddo
  total_volume=0.0_real64;mean_exact=0.0_real64
  do k=1,nz;do j=1,ny;do i=1,nx
    volume=dx(i)*dy(j)*dz(k);total_volume=total_volume+volume
    mean_exact=mean_exact+volume*exact(i,j,k)
  enddo;enddo;enddo
  exact=exact-mean_exact/total_volume
  do k=1,nz;do j=1,ny
    xlo(j,k)=gradx(xf(0));xhi(j,k)=gradx(xf(nx))
  enddo;enddo
  do k=1,nz;do i=1,nx
    ylo(i,k)=grady(yf(0));yhi(i,k)=grady(yf(ny))
  enddo;enddo
  do j=1,ny;do i=1,nx
    zlo(i,j)=gradz(zf(0));zhi(i,j)=gradz(zf(nz))
  enddo;enddo
  do k=1,nz;do j=1,ny;do i=1,nx
    if(i==1)then;gw=xlo(j,k);else;gw=(exact(i,j,k)-exact(i-1,j,k))/(xc(i)-xc(i-1));endif
    if(i==nx)then;ge=xhi(j,k);else;ge=(exact(i+1,j,k)-exact(i,j,k))/(xc(i+1)-xc(i));endif
    if(j==1)then;gs=ylo(i,k);else;gs=(exact(i,j,k)-exact(i,j-1,k))/(yc(j)-yc(j-1));endif
    if(j==ny)then;gn=yhi(i,k);else;gn=(exact(i,j+1,k)-exact(i,j,k))/(yc(j+1)-yc(j));endif
    if(k==1)then;gb=zlo(i,j);else;gb=(exact(i,j,k)-exact(i,j,k-1))/(zc(k)-zc(k-1));endif
    if(k==nz)then;gt=zhi(i,j);else;gt=(exact(i,j,k+1)-exact(i,j,k))/(zc(k+1)-zc(k));endif
    rhs(i,j,k)=(ge-gw)/dx(i)+(gn-gs)/dy(j)+(gt-gb)/dz(k)
  enddo;enddo;enddo

  txlo=xlo;txhi=xhi;tylo=ylo;tyhi=yhi;tzlo=zlo;tzhi=zhi
  call enforce_neumann_compatibility_nonuniform_3d(rhs+0.37_real64,txlo,txhi,tylo,tyhi,tzlo,tzhi, &
       dx,dy,dz,nx,ny,nz,correction,imbalance)
  phi=0.0_real64
  call solve_neumann_poisson_nonuniform_3d(phi,rhs,xlo,xhi,ylo,yhi,zlo,zhi,xc,yc,zc,dx,dy,dz, &
       nx,ny,nz,tol,100000,25,iterations,residual)

  u=1.0_real64;v=2.0_real64;w=3.0_real64
  ox=0.5_real64;oy=-0.25_real64;oz=0.75_real64
  ax=0.2_real64;ay=-0.1_real64;az=0.3_real64;fluid=1
  fluid(nx/2,ny/2,nz/2)=0
  bounds=[xc(2)-epsilon(1.0_real64),xc(nx-1)+epsilon(1.0_real64), &
       yc(2)-epsilon(1.0_real64),yc(ny-1)+epsilon(1.0_real64), &
       zc(2)-epsilon(1.0_real64),zc(nz-1)+epsilon(1.0_real64)]
  call compute_fpm_terms_nonuniform_3d(phi,u,v,w,ox,oy,oz,ax,ay,az,fluid,xc,yc,zc,dx,dy,dz, &
       xlo,xhi,ylo,yhi,zlo,zhi,bounds,gx,gy,gz,vf,af,totals,nx,ny,nz)

  solution_error=maxval(abs(phi-exact));gradient_error=0.0_real64;field_error=0.0_real64
  expected_totals=0.0_real64
  do k=1,nz;do j=1,ny;do i=1,nx
    if(i==1)then;gw=gradx(xf(0));else;gw=gradx(0.5_real64*(xc(i-1)+xc(i)));endif
    if(i==nx)then;ge=gradx(xf(nx));else;ge=gradx(0.5_real64*(xc(i)+xc(i+1)));endif
    if(j==1)then;gs=grady(yf(0));else;gs=grady(0.5_real64*(yc(j-1)+yc(j)));endif
    if(j==ny)then;gn=grady(yf(ny));else;gn=grady(0.5_real64*(yc(j)+yc(j+1)));endif
    if(k==1)then;gb=gradz(zf(0));else;gb=gradz(0.5_real64*(zc(k-1)+zc(k)));endif
    if(k==nz)then;gt=gradz(zf(nz));else;gt=gradz(0.5_real64*(zc(k)+zc(k+1)));endif
    egx=0.5_real64*(gw+ge);egy=0.5_real64*(gs+gn);egz=0.5_real64*(gb+gt)
    cx=3.0_real64*(-0.25_real64)-2.0_real64*0.75_real64
    cy=1.0_real64*0.75_real64-3.0_real64*0.5_real64
    cz=2.0_real64*0.5_real64-1.0_real64*(-0.25_real64)
    ev=-real(fluid(i,j,k),real64)*(cx*egx+cy*egy+cz*egz)
    ea=-real(fluid(i,j,k),real64)*(0.2_real64*egx-0.1_real64*egy+0.3_real64*egz)
    gradient_error=max(gradient_error,abs(gx(i,j,k)-egx),abs(gy(i,j,k)-egy),abs(gz(i,j,k)-egz))
    field_error=max(field_error,abs(vf(i,j,k)-ev),abs(af(i,j,k)-ea))
    inside=merge(1.0_real64,0.0_real64,xc(i)>bounds(1).and.xc(i)<bounds(2).and. &
         yc(j)>bounds(3).and.yc(j)<bounds(4).and.zc(k)>bounds(5).and.zc(k)<bounds(6))
    volume=inside*dx(i)*dy(j)*dz(k)
    expected_totals=expected_totals+volume*[ev,ea]
  enddo;enddo;enddo
  integral_error=maxval(abs(totals-expected_totals))
  write(output_unit,'(a,es14.6)')'compatibility correction: ',correction
  write(output_unit,'(a,es14.6)')'compatibility imbalance: ',imbalance
  write(output_unit,'(a,i0,a,es14.6)')'iterations: ',iterations,', residual: ',residual
  write(output_unit,'(a,es14.6)')'solution error (Linf): ',solution_error
  write(output_unit,'(a,es14.6)')'gradient error (Linf): ',gradient_error
  write(output_unit,'(a,es14.6)')'FPM field error (Linf): ',field_error
  write(output_unit,'(a,es14.6)')'integral error (Linf): ',integral_error
  write(output_unit,'(a,2(es14.6,1x))')'integrated terms: ',totals
  passed=iterations<100000.and.residual<=tol.and.abs(correction-0.37_real64/6.0_real64)<1e-12
  passed=passed.and.abs(imbalance)<1e-12.and.solution_error<1e-8.and.gradient_error<1e-8
  passed=passed.and.field_error<1e-8.and.integral_error<1e-10
  if(.not.passed)then;write(error_unit,'(a)')'FAIL: nonuniform 3-D Laplace/FPM test';error stop 1;endif
  write(output_unit,'(a)')'PASS: nonuniform 3-D Laplace/FPM test'
contains
  subroutine make_axis(face,center,width,n,amplitude)
    integer,intent(in)::n
    real(real64),intent(out)::face(0:n),center(n),width(n)
    real(real64),intent(in)::amplitude
    real(real64)::weight(n),scale
    integer::m
    do m=1,n;weight(m)=1.0_real64+amplitude*sin(2.0_real64*acos(-1.0_real64)*(m-0.5_real64)/n);enddo
    scale=1.0_real64/sum(weight);face(0)=0.0_real64
    do m=1,n;width(m)=scale*weight(m);face(m)=face(m-1)+width(m);center(m)=0.5_real64*(face(m-1)+face(m));enddo
  end subroutine make_axis
  pure function field(x,y,z) result(value)
    real(real64),intent(in)::x,y,z;real(real64)::value
    value=0.2_real64*x*x+0.3_real64*y*y-0.1_real64*z*z+0.4_real64*x-0.2_real64*y+0.3_real64*z
  end function field
  pure function gradx(x) result(value)
    real(real64),intent(in)::x;real(real64)::value;value=0.4_real64*x+0.4_real64
  end function gradx
  pure function grady(y) result(value)
    real(real64),intent(in)::y;real(real64)::value;value=0.6_real64*y-0.2_real64
  end function grady
  pure function gradz(z) result(value)
    real(real64),intent(in)::z;real(real64)::value;value=-0.2_real64*z+0.3_real64
  end function gradz
end program test_fpm_3d_nonuniform
