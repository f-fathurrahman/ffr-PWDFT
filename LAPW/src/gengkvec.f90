! !DESCRIPTION:
!   Generates a set of ${\bf G+k}$-vectors for the input $k$-point with length
!   less than {\tt gkmax}.
SUBROUTINE gengkvec(ngvec,ivg,vgc,vkl,vkc,gkmax,ngkmax,ngk,igkig,vgkl,vgkc,gkc)
! !INPUT/OUTPUT PARAMETERS:
!   ngvec  : number of G-vectors (in,integer)
!   ivg    : G-vector INTEGER coordinates (in,integer(3,ngvec))
!   vgc    : G-vectors in Cartesian coordinates (in,real(3,ngvec))
!   vkl    : k-point vector in lattice coordinates (in,real(3))
!   vkc    : k-point vector in Cartesian coordinates (in,real(3))
!   gkmax  : G+k-vector cut-off (in,real)
!   ngkmax : maximum number of G+k-vectors (in,integer)
!   ngk    : number of G+k-vectors RETURN ed (out,integer)
!   igkig  : index from G+k-vectors to G-vectors (out,integer(ngkmax))
!   vgkl   : G+k-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgkc   : G+k-vectors in Cartesian coordinates (out,real(3,ngkmax))
!   gkc    : length of G+k-vectors (out,real(ngkmax))
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngvec,ivg(3,ngvec)
  REAL(8), intent(in) :: vgc(3,ngvec)
  REAL(8), intent(in) :: vkl(3),vkc(3)
  REAL(8), intent(in) :: gkmax
  INTEGER, intent(in) :: ngkmax
  INTEGER, intent(out) :: ngk,igkig(ngkmax)
  REAL(8), intent(out) :: vgkl(3,ngkmax),vgkc(3,ngkmax),gkc(ngkmax)
  ! local variables
  INTEGER :: ig
  REAL(8) :: v1,v2,v3,t0,t1
  
  t0=gkmax**2
  ngk=0
  DO ig=1,ngvec
    v1 = vgc(1,ig) + vkc(1)
    v2 = vgc(2,ig) + vkc(2)
    v3 = vgc(3,ig) + vkc(3)
    t1 = v1**2 + v2**2 + v3**2
    IF(t1 < t0) THEN 
      ngk = ngk + 1
      ! index to G-vector
      igkig(ngk) = ig
      ! G+k-vector in lattice coordinates
      vgkl(:,ngk) = dble(ivg(:,ig))+vkl(:)
      ! G+k-vector in Cartesian coordinates
      vgkc(1,ngk) = v1
      vgkc(2,ngk) = v2
      vgkc(3,ngk) = v3
      ! length of G+k-vector
      gkc(ngk) = sqrt(t1)
      IF(ngk==ngkmax) exit
    ENDIF 
  ENDDO 
  RETURN 
END SUBROUTINE 
