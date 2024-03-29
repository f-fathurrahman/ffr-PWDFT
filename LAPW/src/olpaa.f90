SUBROUTINE olpaa(tor,ias,ngp,apwalm,ld,o)
  USE m_atoms, ONLY: idxis
  USE m_gkvectors, ONLY: ngkmax
  USE m_apwlo, ONLY: apwordmax, apword, lmoapw
  USE m_muffin_tins, ONLY: lmaxapw, lmmaxapw
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: tor
  INTEGER, intent(in) :: ias,ngp
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(inout) :: o(*)
  ! local variables
  INTEGER :: is,lmo,io
  INTEGER :: l,m,lm,i
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: a(:,:)

  is = idxis(ias)
  lmo = lmoapw(is)
  ALLOCATE( a(lmo,ngp) )
  i = 0
  lm = 0
  DO l=0,lmaxapw
    DO m=-l,l
      lm=lm+1
      DO io=1,apword(l,is)
        i = i+1
        a(i,1:ngp) = apwalm(1:ngp,io,lm)
      ENDDO 
    ENDDO 
  ENDDO 
  CALL zmctmu(tor,lmo,ngp,a,a,ld,o)
  DEALLOCATE(a)
  RETURN 
END SUBROUTINE 

