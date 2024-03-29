! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   ld    : leading dimension of o (in,integer)
!   o     : overlap matrix (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
SUBROUTINE olpistl(ngp,igpig,ld,o)
  USE m_gkvectors, ONLY: ngkmax
  USE m_gvectors, ONLY: cfunig, ivgig, ivg
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp,igpig(ngkmax)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(inout) :: o(*)
  ! local variables
  INTEGER :: i1,i2,i3,j1,j2,j3
  INTEGER :: ig,i,j,k
  DO j=1,ngp
    k=(j-1)*ld
    ig=igpig(j)
    j1=ivg(1,ig)
    j2=ivg(2,ig)
    j3=ivg(3,ig)
    DO i=1,j
      k=k+1
      ig=igpig(i)
      i1=ivg(1,ig)-j1
      i2=ivg(2,ig)-j2
      i3=ivg(3,ig)-j3
      o(k) = o(k) + cfunig(ivgig(i1,i2,i3))
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 
