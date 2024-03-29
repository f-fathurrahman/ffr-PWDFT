! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   zfmt1 : input complex muffin-tin function in spherical coordinates
!           (in,complex(*))
!   zfmt2 : output complex muffin-tin function in spherical harmonics
!           (out,complex(*))
SUBROUTINE zfsht(nr,nri,zfmt1,zfmt2)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo
  USE m_sht, ONLY: zfshti, zfshto
  USE m_constants, ONLY: zone, zzero
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri
  COMPLEX(8), intent(in) :: zfmt1(*)
  COMPLEX(8), intent(out) :: zfmt2(*)
  ! local variables
  INTEGER nro,i
  ! transform the inner part of the muffin-tin
  CALL zgemm('N','N',lmmaxi,nri,lmmaxi,zone,zfshti,lmmaxi,zfmt1,lmmaxi,zzero, &
   zfmt2,lmmaxi)
  ! transform the outer part of the muffin-tin
  nro=nr-nri
  i=lmmaxi*nri+1
  CALL zgemm('N','N',lmmaxo,nro,lmmaxo,zone,zfshto,lmmaxo,zfmt1(i),lmmaxo,zzero, &
   zfmt2(i),lmmaxo)
  RETURN 
END SUBROUTINE 
!EOC

