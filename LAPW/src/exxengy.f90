SUBROUTINE exxengy()
  USE modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,ist,jst,is,ia
  INTEGER :: nrc,nrci,npc
  INTEGER :: m1,m2
  COMPLEX(8) :: z1
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: wfcr1(:,:),wfcr2(:,:)
  COMPLEX(8), ALLOCATABLE :: zrhomt(:),zvclmt(:),zfmt(:)
  ! external functions
  COMPLEX(8) zfmtinp
  external zfmtinp

  ALLOCATE(wfcr1(npcmtmax,2),wfcr2(npcmtmax,2))
  ALLOCATE(zrhomt(npcmtmax),zvclmt(npcmtmax),zfmt(npcmtmax))
  ! zero the exchange energy
  engyx=0.d0

  !--------------------------------------------------!
  !     val-val-val and val-cr-val contributions     !
  !--------------------------------------------------!
  DO ik=1,nkpt
    WRITE(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
    CALL exxengyk(ik)
  ENDDO 

  !-----------------------------------!
  !    core-core-core contribution    !
  !-----------------------------------!
  ! begin loops over atoms and species
  DO is=1,nspecies
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    DO ia=1,natoms(is)
      DO jst=1,nstsp(is)
        IF(spcore(jst,is)) THEN 
          DO m2=-ksp(jst,is),ksp(jst,is)-1
            ! generate the core wavefunction in spherical coordinates (pass in m-1/2)
            CALL wavefcr(.false.,lradstp,is,ia,jst,m2,npcmtmax,wfcr2)
            DO ist=1,nstsp(is)
              IF(spcore(ist,is)) THEN 
                DO m1=-ksp(ist,is),ksp(ist,is)-1
                  CALL wavefcr(.false.,lradstp,is,ia,ist,m1,npcmtmax,wfcr1)
                  ! calculate the complex overlap density
                  CALL zrho2(npc,wfcr1,wfcr1(:,2),wfcr2,wfcr2(:,2),zfmt)
                  CALL zfsht(nrc,nrci,zfmt,zrhomt)
                  ! calculate the Coulomb potential
                  CALL zpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
                   zrhomt,zvclmt)
                  z1=zfmtinp(nrc,nrci,wrcmt(:,is),zrhomt,zvclmt)
                  engyx=engyx-0.5d0*dble(z1)
                ENDDO 
              ENDIF 
            ENDDO ! ist
          ENDDO  ! jst
        ENDIF 
      ENDDO 
    ! end loops over atoms and species
    ENDDO 
  ENDDO 
  DEALLOCATE(wfcr1,wfcr2,zrhomt,zvclmt,zfmt)
  RETURN 

CONTAINS 

SUBROUTINE zrho2(n,x1,x2,y1,y2,z)
  IMPLICIT NONE 
  INTEGER, intent(in) :: n
  COMPLEX(8), intent(in) :: x1(n),x2(n),y1(n),y2(n)
  COMPLEX(8), intent(out) :: z(n)
  z(:)=conjg(x1(:))*y1(:)+conjg(x2(:))*y2(:)
  RETURN 
END SUBROUTINE 

END SUBROUTINE 

