! !DESCRIPTION:
!   Calculates and stores the spherical Bessel functions
!   $j_l(|{\bf G}+{\bf p}|{\bf R}_{\rm MT})$ for all input ${\bf G}+{\bf p}$
!   vectors and the muffin-tin radii ${\bf R}_{\rm MT}$ of every atomic species.
SUBROUTINE genjlgprmt(lmax,ngp,gpc,ld,jlgprmt)
  USE m_atoms, ONLY: nspecies
  USE m_muffin_tins, ONLY: rmt
! !INPUT/OUTPUT PARAMETERS:
!   lmax    : angular momentum cut-off (in,integer)
!   ngp     : number of G+p-vectors (in,integer)
!   gpc     : length of G+p-vectors (in,real(ngkmax))
!   ld      : leading dimension (in,integer)
!   jlgprmt : spherical Bessel functions (out,real(0:lmax,ld,nspecies))
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lmax,ngp
  REAL(8), intent(in) :: gpc(ngp)
  INTEGER, intent(in) :: ld
  REAL(8), intent(out) :: jlgprmt(0:lmax,ld,nspecies)
  ! local variables
  INTEGER is,ig
  REAL(8) r,t1
  DO is=1,nspecies
    r=rmt(is)
    DO ig=1,ngp
      t1=gpc(ig)*r
      CALL sbessel(lmax,t1,jlgprmt(:,ig,is))
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 
