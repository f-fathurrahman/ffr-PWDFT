
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potefield
use modmain
implicit none
! local variables
integer is,ia,ias
integer nr,nri,ir
integer i,i1,i2,i3
real(8) e,tp(2),r,t1
real(8) v0,e00,elm(-1:1)
real(8) v1(3),v2(3)
! constant added to potential so that it is zero at the unit cell center
v1(:)=0.5d0*(avec(:,1)+avec(:,2)+avec(:,3))
v0=dot_product(efieldc(:),v1(:))
! determine the electric field vector in spherical coordinates
call sphcrd(efieldc,e,tp)
! coefficients for real spherical harmonics R_1-1, R_10 and R_11
t1=e*sqrt(fourpi/3.d0)
elm(-1)=t1*sin(tp(1))*sin(tp(2))
elm(0)=-t1*cos(tp(1))
elm(1)=t1*sin(tp(1))*cos(tp(2))
! muffin-tin potential
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! coefficient for R_00
    e00=v0-dot_product(efieldc(:),atposc(:,ia,is))
    e00=e00/y00
    i=1
    do ir=1,nr
      r=rsp(ir,is)
      vclmt(i,ias)=vclmt(i,ias)+e00
      vclmt(i+1,ias)=vclmt(i+1,ias)+elm(-1)*r
      vclmt(i+2,ias)=vclmt(i+2,ias)+elm(0)*r
      vclmt(i+3,ias)=vclmt(i+3,ias)+elm(1)*r
      if (ir.le.nri) then
        i=i+lmmaxi
      else
        i=i+lmmaxo
      end if
    end do
  end do
end do
! interstitial potential
ir=0
do i3=0,ngridg(3)-1
  v1(3)=dble(i3)/dble(ngridg(3))
  do i2=0,ngridg(2)-1
    v1(2)=dble(i2)/dble(ngridg(2))
    do i1=0,ngridg(1)-1
      v1(1)=dble(i1)/dble(ngridg(1))
      ir=ir+1
      call r3mv(avec,v1,v2)
      vclir(ir)=vclir(ir)+v0-dot_product(efieldc(:),v2(:))
    end do
  end do
end do
return
end subroutine

