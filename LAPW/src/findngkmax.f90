! !DESCRIPTION:
!   Determines the largest number of ${\bf G+k}$-vectors with length less than
!   {\tt gkmax} over all the $k$-points. This variable is used for allocating
!   arrays.
SUBROUTINE findngkmax(nkpt,vkc,nspnfv,vqcss,ngvec,vgc,gkmax,ngkmax)
! !INPUT/OUTPUT PARAMETERS:
!   nkpt   : number of k-points (in,integer)
!   vkc    : k-point vectors in Cartesian coordinates (in,real(3,nkpt))
!   nspnfv : number of first-variational spin components: 1 normal case, 2 for
!            spin-spiral case (in,integer)
!   vqcss  : spin-spiral q-vector, not referenced if nspnfv=1 (in,integer)
!   ngvec  : number of G-vectors (in,integer)
!   vgc    : G-vectors in Cartesian coordinates (in,real(3,ngvec))
!   gkmax  : maximum allowed |G+k| (in,real)
!   ngkmax : maximum number of G+k-vectors over all k-points (out,integer)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nkpt
  REAL(8), intent(in) :: vkc(3,nkpt)
  INTEGER, intent(in) :: nspnfv
  REAL(8), intent(in) :: vqcss(3)
  INTEGER, intent(in) :: ngvec
  REAL(8), intent(in) :: vgc(3,ngvec)
  REAL(8), intent(in) :: gkmax
  INTEGER, intent(out) :: ngkmax
  ! local variables
  INTEGER ispn,ik,n,ig
  REAL(8) v1(3),v2(3),t0,t1
  t0=gkmax**2+1.d-6
  ngkmax=0
  DO ispn=1,nspnfv
    DO ik=1,nkpt
      IF(nspnfv==2) THEN 
        ! spin-spiral case
        IF(ispn==1) THEN 
          v1(:) = vkc(:,ik) + 0.5d0*vqcss(:)
        ELSE 
          v1(:) = vkc(:,ik) - 0.5d0*vqcss(:)
        ENDIF 
      ELSE 
        v1(:) = vkc(:,ik)
      ENDIF 
      n=0
      DO ig = 1,ngvec
        v2(:) = vgc(:,ig) + v1(:)
        t1 = v2(1)**2 + v2(2)**2 + v2(3)**2
        IF(t1 < t0) n=n+1
      ENDDO 
      ngkmax = max(ngkmax,n)
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 
