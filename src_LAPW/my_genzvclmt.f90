SUBROUTINE my_genzvclmt(nr,nri,ld1,rl,wpr,ld2,zrhomt,zvclmt)
! output: zvclmt
  USE m_atoms, ONLY: nspecies, natmtot, idxis
  USE m_muffin_tins, ONLY: lmaxo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies)
  INTEGER, intent(in) :: ld1
  REAL(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
  REAL(8), intent(in) :: wpr(4,ld1,nspecies)
  INTEGER, intent(in) :: ld2
  COMPLEX(8), intent(in) :: zrhomt(ld2,natmtot)
  COMPLEX(8), intent(out) :: zvclmt(ld2,natmtot)
  ! local variables
  INTEGER :: is,ias

  write(*,*)
  write(*,*) 'In my_genzvclmt: ld1 = ', ld1
  write(*,*) 'In my_genzvclmt: ld2 = ', ld2  

  ! Loop over all atoms, results accumulated in zvclmt
  DO ias = 1,natmtot
    is = idxis(ias)
    !
    zvclmt(:,ias) = cmplx(0.d0, 0.d0, kind=8)
    write(*,'(1x,A,2ES18.10)') 'in my_genzvclmt: sum(zrhomt[ias]) before = ', sum(zrhomt(1:ld2,ias))
    write(*,'(1x,A,2ES18.10)') 'in my_genzvclmt: sum(zvclmt[ias]) before = ', sum(zvclmt(1:ld2,ias))
    !
    CALL my_zpotclmt( nr(is), nri(is), ld1, rl(:,:,is), wpr(:,:,is), &
      zrhomt(:,ias), zvclmt(:,ias) )
    !
    write(*,'(1x,A,2ES18.10)') 'in my_genzvclmt: sum(zrhomt[ias]) after = ', sum(zrhomt(1:ld2,ias))
    write(*,'(1x,A,2ES18.10)') 'in my_genzvclmt: sum(zvclmt[ias]) after = ', sum(zvclmt(1:ld2,ias))
  ENDDO 

  RETURN 
END SUBROUTINE 
