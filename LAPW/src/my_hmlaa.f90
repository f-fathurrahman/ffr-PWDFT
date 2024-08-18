! !INPUT/OUTPUT PARAMETERS:
!   thr    : .true. if the matrix h is real valued (in,logical)
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   ld     : leading dimension of h (in,integer)
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
SUBROUTINE my_hmlaa(thr,ias,ngp,apwalm,ld,h)
  USE m_atoms, ONLY: idxis
  USE m_gkvectors, ONLY: ngkmax
  USE m_apwlo, ONLY: apwordmax, apwfr, apwdfr, apword, lmoapw
  USE m_muffin_tins, ONLY: nrmt, idxlm, rmt, lmaxo, lmaxapw, lmmaxapw
  USE m_hamiltonian, ONLY: haa, gntyry
  IMPLICIT NONE 
  ! arguments
  LOGICAL, INTENT(in) :: thr
  INTEGER, INTENT(in) :: ias,ngp
  COMPLEX(8), INTENT(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
  INTEGER, INTENT(in) :: ld
  COMPLEX(8), INTENT(inout) :: h(*)
  ! local variables
  INTEGER :: is,lmo,io,jo,i
  INTEGER :: l1,l2,l3,m1,m2,m3
  INTEGER :: lm1,lm2,lm3
  REAL(8) :: t0
  COMPLEX(8) :: z1
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: a(:,:),b(:,:)

  is = idxis(ias)
  lmo = lmoapw(is)
  !
  ALLOCATE(a(lmo,ngp),b(lmo,ngp))
  !
  t0 = 0.5d0*rmt(is)**2
  !
  i = 0
  lm1 = 0
  DO l1 = 0,lmaxapw
    DO m1 = -l1,l1
      lm1 = lm1 + 1
      DO io = 1,apword(l1,is)
        i = i + 1
        b(i,:) = 0.d0
        lm3 = 0
        DO l3 = 0,lmaxapw
          DO m3 = -l3,l3
            lm3 = lm3 + 1
            DO jo = 1,apword(l3,is)
              z1 = 0.d0 ! complex
              DO l2 = 0,lmaxo
                IF( mod(l1+l2+l3,2) == 0 ) THEN 
                  DO m2=-l2,l2
                    lm2 = idxlm(l2,m2)
                    z1 = z1 + gntyry(lm2,lm3,lm1)*haa(lm2,jo,l3,io,l1,ias)
                  ENDDO 
                ENDIF 
              ENDDO 
              IF( abs(dble(z1)) + abs(aimag(z1)) > 1.d-14 ) THEN 
                !write(*,*) 'z1 = ', z1
                call zaxpy(ngp, z1, apwalm(:,jo,lm3),1, b(i,1),lmo )
              ENDIF 
            ENDDO 
          ENDDO 
        ENDDO 
        ! kinetic surface contribution
        DO jo = 1,apword(l1,is)
          z1 = t0*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
          !write(*,*) 'In kinetic surface contribution: z1 = ', z1
          CALL zaxpy(ngp,z1,apwalm(:,jo,lm1),1,b(i,1),lmo)
        ENDDO 
        a(i,1:ngp) = apwalm(1:ngp,io,lm1)
      ENDDO 
    ENDDO 
  ENDDO 
  !write(*,*) 'Before zmctmu: sum(a) = ', sum(a)
  !write(*,*) 'Before zmctmu: sum(b) = ', sum(b)
  CALL zmctmu(thr, lmo, ngp, a, b, ld, h)
  DEALLOCATE(a,b)
  RETURN 
END SUBROUTINE 
