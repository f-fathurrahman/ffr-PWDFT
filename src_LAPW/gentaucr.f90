SUBROUTINE gentaucr()
  USE modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER ist,ispn,jspn
  INTEGER is,ia,ias,nthd
  INTEGER nr,nri,np,i,m
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)
  COMPLEX(8), ALLOCATABLE :: wfcr(:,:),gzfmt(:,:),zfmt(:)

  taucr(:,:,:)=0.d0
  ALLOCATE(wfcr(npmtmax,2),gzfmt(npmtmax,3),zfmt(npmtmax))
  !
  DO ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    np=npmt(is)
    !
    DO ist=1,nstsp(is)
      IF(spcore(ist,is)) THEN 
        DO m=-ksp(ist,is),ksp(ist,is)-1
          ! generate the core wavefunction in spherical harmonics (pass in m-1/2)
          CALL wavefcr(.true.,1,is,ia,ist,m,npmtmax,wfcr)
          DO ispn=1,2
            IF(spinpol) THEN 
              jspn=ispn
            ELSE 
              jspn=1
            ENDIF 
            ! compute the gradient
            CALL gradzfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),wfcr(:,ispn), &
             npmtmax,gzfmt)
            DO i=1,3
              ! convert gradient to spherical coordinates
              CALL zbsht(nr,nri,gzfmt(:,i),zfmt)
              ! add to total in muffin-tin
              taucr(1:np,ias,jspn)=taucr(1:np,ias,jspn) + &
                  0.5d0*(dble(zfmt(1:np))**2+aimag(zfmt(1:np))**2)
            ENDDO 
          ENDDO 
        ENDDO 
      ENDIF 
    ENDDO 
  ENDDO 
  DEALLOCATE(wfcr,gzfmt,zfmt)

  ! convert core tau to spherical harmonics
  ALLOCATE(rfmt(npmtmax))
  DO ispn=1,nspinor
    DO ias=1,natmtot
      is=idxis(ias)
      CALL dcopy(npmt(is),taucr(:,ias,ispn),1,rfmt,1)
      CALL rfsht(nrmt(is),nrmti(is),rfmt,taucr(:,ias,ispn))
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt)
  RETURN 

END SUBROUTINE 

