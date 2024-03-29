! !DESCRIPTION:
!   Generates the atomic structure factors for a set of ${\bf G+p}$-vectors:
!   $$ S_{\alpha}({\bf G+p})=\exp(i({\bf G+p})\cdot{\bf r}_{\alpha}), $$
!   where ${\bf r}_{\alpha}$ is the position of atom $\alpha$.
SUBROUTINE gensfacgp(ngp,vgpc,ld,sfacgp)
  ! !USES:
  USE m_atoms, ONLY: atposc, idxia, idxis, natmtot
  ! !INPUT/OUTPUT PARAMETERS:
  !   ngp    : number of G+p-vectors (in,integer)
  !   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,*))
  !   ld     : leading dimension (in,integer)
  !   sfacgp : structure factors of G+p-vectors (out,complex(ld,natmtot))
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp
  REAL(8), intent(in) :: vgpc(3,ngp)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(out) :: sfacgp(ld,natmtot)
  ! local variables
  INTEGER :: is,ia,ias
  INTEGER :: igp
  REAL(8) :: v1,v2,v3,t1
  
  DO ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    v1=atposc(1,ia,is); v2=atposc(2,ia,is); v3=atposc(3,ia,is)
    DO igp=1,ngp
      t1=vgpc(1,igp)*v1+vgpc(2,igp)*v2+vgpc(3,igp)*v3
      sfacgp(igp,ias)=cmplx(cos(t1),sin(t1),8)
    ENDDO 
  ENDDO 
  RETURN 

END SUBROUTINE 
