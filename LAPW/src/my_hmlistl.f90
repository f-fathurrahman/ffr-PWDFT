! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   ld    : leading dimension of h (in,integer)
!   h     : Hamiltonian matrix (inout,complex(*))
SUBROUTINE my_hmlistl(ngp, igpig, vgpc, ld, h)
  USE m_gkvectors, ONLY: ngkmax
  USE m_gvectors, ONLY: ivgig, cfunig, ivg
  USE m_density_pot_xc, ONLY: vsig 
  IMPLICIT NONE
  ! arguments
  INTEGER, INTENT(in) :: ngp,igpig(ngkmax)
  REAL(8), INTENT(in) :: vgpc(3,ngkmax)
  INTEGER, INTENT(in) :: ld
  COMPLEX(8), INTENT(inout) :: h(*)
  ! local variables
  INTEGER :: i1,i2,i3,j1,j2,j3
  INTEGER :: ig,i,j,k
  REAL(8) :: v1,v2,v3,t1
  
  DO j = 1,ngp
    k = (j-1)*ld
    ig = igpig(j)
    j1 = ivg(1,ig)
    j2 = ivg(2,ig)
    j3 = ivg(3,ig)
    v1 = 0.5d0*vgpc(1,j)
    v2 = 0.5d0*vgpc(2,j)
    v3 = 0.5d0*vgpc(3,j)
    DO i = 1,j
      k = k + 1
      ig = igpig(i)
      i1 = ivg(1,ig) - j1
      i2 = ivg(2,ig) - j2
      i3 = ivg(3,ig) - j3
      ig = ivgig(i1,i2,i3)
      t1 = vgpc(1,i)*v1 + vgpc(2,i)*v2 + vgpc(3,i)*v3
      h(k) = h(k) + vsig(ig) + t1*cfunig(ig)
    ENDDO 
  ENDDO 
  RETURN 
  
END SUBROUTINE 

