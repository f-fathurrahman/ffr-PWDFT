
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfplot(np,vpl,rfmt,rfir,fp)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: np
real(8), intent(in) :: vpl(3,np)
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
real(8), intent(out) :: fp(np)
! local variables
integer ias,is,ip,nthd
! allocatable arrays
real(8), allocatable :: rfmt1(:,:,:)
complex(8), allocatable :: zfft(:)
! unpack the muffin-tin function
allocate(rfmt1(lmmaxo,nrmtmax,natmtot))
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call rfmtpack(.false.,nrmt(is),nrmti(is),rfmt(:,ias),rfmt1(:,:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! Fourier transform rfir to G-space
allocate(zfft(ngtot))
zfft(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft)
! begin loop over all points
call holdthd(np,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ip=1,np
  call rfip(ip)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(rfmt1,zfft)
return

contains

subroutine rfip(ip)
implicit none
! arguments
integer, intent(in) :: ip
! local variables
integer is,ia,ias,nr,nri
integer ir0,ir,lmax,l,m,lm
integer ig,ifg,i1,i2,i3,i,j
real(8) rmt2,r,sum,ya(4),t1
real(8) v1(3),v2(3),v3(3),v4(3),v5(3)
! automatic arrays
real(8) rlm(lmmaxo)
v2(:)=vpl(:,ip)
call r3frac(epslat,v2)
! convert point to Cartesian coordinates
call r3mv(avec,v2,v1)
! check if point is in a muffin-tin
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    v2(:)=v1(:)-atposc(:,ia,is)
    do i1=-1,1
      v3(:)=v2(:)+dble(i1)*avec(:,1)
      do i2=-1,1
        v4(:)=v3(:)+dble(i2)*avec(:,2)
        do i3=-1,1
          v5(:)=v4(:)+dble(i3)*avec(:,3)
          t1=v5(1)**2+v5(2)**2+v5(3)**2
          if (t1.lt.rmt2) then
            r=sqrt(t1)
            call genrlmv(lmaxo,v5,rlm)
            do ir=1,nr
              if (rsp(ir,is).ge.r) then
                if (ir.le.3) then
                  ir0=1
                else if (ir.gt.nr-2) then
                  ir0=nr-3
                else
                  ir0=ir-2
                end if
                r=max(r,rsp(1,is))
                if (ir0.le.nri) then
                  lmax=lmaxi
                else
                  lmax=lmaxo
                end if
                sum=0.d0
                lm=0
                do l=0,lmax
                  do m=-l,l
                    lm=lm+1
                    do j=1,4
                      i=ir0+j-1
                      ya(j)=rfmt1(lm,i,ias)
                    end do
                    t1=poly4(rsp(ir0,is),ya,r)
                    sum=sum+t1*rlm(lm)
                  end do
                end do
                goto 10
              end if
            end do
          end if
        end do
      end do
    end do
  end do
end do
! otherwise use direct Fourier transform of interstitial function
sum=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  t1=vgc(1,ig)*v1(1)+vgc(2,ig)*v1(2)+vgc(3,ig)*v1(3)
  sum=sum+dble(zfft(ifg)*cmplx(cos(t1),sin(t1),8))
end do
10 continue
fp(ip)=sum
return
end subroutine

real(8) function poly4(xa,ya,x)
implicit none
! arguments
real(8), intent(in) :: xa(4),ya(4),x
! local variables
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
! evaluate the polynomial coefficients
x0=xa(1)
x1=xa(2)-x0; x2=xa(3)-x0; x3=xa(4)-x0
t4=x1-x2; t5=x1-x3; t6=x2-x3
y0=ya(1)
y1=ya(2)-y0; y2=ya(3)-y0; y3=ya(4)-y0
t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
t0=1.d0/(x2*t3*t4*t5*t6)
t3=t3*y2
c3=t1*t4+t2*t6-t3*t5
t4=x1**2; t5=x2**2; t6=x3**2
c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
t1=x-x0
! evaluate the polynomial
poly4=y0+t0*t1*(c1+t1*(c2+c3*t1))
return
end function

end subroutine

