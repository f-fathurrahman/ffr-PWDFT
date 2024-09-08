
! Copyright (C) 2002-2012 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengkvec
! !INTERFACE:
subroutine gengkvec(ngvec,ivg,vgc,vkl,vkc,gkmax,ngkmax,ngk,igkig,vgkl,vgkc,gkc)
! !INPUT/OUTPUT PARAMETERS:
!   ngvec  : number of G-vectors (in,integer)
!   ivg    : G-vector integer coordinates (in,integer(3,ngvec))
!   vgc    : G-vectors in Cartesian coordinates (in,real(3,ngvec))
!   vkl    : k-point vector in lattice coordinates (in,real(3))
!   vkc    : k-point vector in Cartesian coordinates (in,real(3))
!   gkmax  : G+k-vector cut-off (in,real)
!   ngkmax : maximum number of G+k-vectors (in,integer)
!   ngk    : number of G+k-vectors returned (out,integer)
!   igkig  : index from G+k-vectors to G-vectors (out,integer(ngkmax))
!   vgkl   : G+k-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgkc   : G+k-vectors in Cartesian coordinates (out,real(3,ngkmax))
!   gkc    : length of G+k-vectors (out,real(ngkmax))
! !DESCRIPTION:
!   Generates a set of ${\bf G+k}$-vectors for the input $k$-point with length
!   less than {\tt gkmax}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed spherical coordinate generation, May 2010 (JKD)
!   Removed modmain and added arguments, September 2012 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngvec,ivg(3,ngvec)
real(8), intent(in) :: vgc(3,ngvec)
real(8), intent(in) :: vkl(3),vkc(3)
real(8), intent(in) :: gkmax
integer, intent(in) :: ngkmax
integer, intent(out) :: ngk,igkig(ngkmax)
real(8), intent(out) :: vgkl(3,ngkmax),vgkc(3,ngkmax),gkc(ngkmax)
! local variables
integer ig
real(8) v1,v2,v3,t0,t1
t0=gkmax**2
ngk=0
do ig=1,ngvec
  v1=vgc(1,ig)+vkc(1)
  v2=vgc(2,ig)+vkc(2)
  v3=vgc(3,ig)+vkc(3)
  t1=v1**2+v2**2+v3**2
  if (t1.lt.t0) then
    ngk=ngk+1
! index to G-vector
    igkig(ngk)=ig
! G+k-vector in lattice coordinates
    vgkl(:,ngk)=dble(ivg(:,ig))+vkl(:)
! G+k-vector in Cartesian coordinates
    vgkc(1,ngk)=v1
    vgkc(2,ngk)=v2
    vgkc(3,ngk)=v3
! length of G+k-vector
    gkc(ngk)=sqrt(t1)
    if (ngk.eq.ngkmax) exit
  end if
end do
return
end subroutine
!EOC

