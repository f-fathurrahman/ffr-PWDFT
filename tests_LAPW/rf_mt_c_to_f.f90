! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot))
SUBROUTINE rf_mt_c_to_f(rfmt)
  USE m_atoms, ONLY: natmtot, rsp, idxis
  USE m_muffin_tins, ONLY: npmtmax, nrcmtmax, nrmtmax, npmt, rcmt, npcmti, nrcmt, nrcmti, &
                   npmti, lmmaxo, lmmaxi, nrmti, nrmt, lradstp, rlmt
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(inout) :: rfmt(npmtmax,natmtot)
  ! local variables
  INTEGER :: is,ias,lm
  INTEGER :: nr,nri,nro
  INTEGER :: iro,ir,npi,i
  INTEGER :: nrc,nrci,nrco
  INTEGER :: irco,irc,npci
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: fi(:),fo(:),rfmt1(:)

  IF(lradstp == 1) RETURN 
  ALLOCATE(fi(nrcmtmax),fo(nrmtmax),rfmt1(npmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    nro=nr-nri
    iro=nri+1
    npi=npmti(is)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    nrco=nrc-nrci
    irco=nrci+1
    npci=npcmti(is)
    ! interpolate up to lmaxi over entire muffin-tin
    DO lm=1,lmmaxi
      i=lm
      DO irc=1,nrci
        fi(irc)=rfmt(i,ias)
        i=i+lmmaxi
      ENDDO 
      DO irc=irco,nrc
        fi(irc)=rfmt(i,ias)
        i=i+lmmaxo
      ENDDO 
      CALL rf_interp(nrc,rcmt(:,is),fi,nr,rlmt(:,1,is),fo)
      i=lm
      DO ir=1,nri
        rfmt1(i)=fo(ir)
        i=i+lmmaxi
      ENDDO 
      DO ir=iro,nr
        rfmt1(i)=fo(ir)
        i=i+lmmaxo
      ENDDO 
    ENDDO 
    ! interpolate up to lmaxo on outer part of muffin-tin
    DO lm=lmmaxi+1,lmmaxo
      i=npci+lm
      DO irc=irco,nrc
        fi(irc)=rfmt(i,ias)
        i=i+lmmaxo
      ENDDO 
      !write(*,*) 'irco = ', irco, ' nrco = ', nrco
      !write(*,*) 'iro = ', iro, ' nro = ', nro
      CALL rf_interp(nrco, rcmt(irco,is), fi(irco),nro,rsp(iro,is),fo(iro))
      !stop 'ffr 65'
      i=npi+lm
      DO ir=iro,nr
        rfmt1(i)=fo(ir)
        i=i+lmmaxo
      ENDDO 
    ENDDO 
    CALL dcopy(npmt(is),rfmt1,1,rfmt(:,ias),1)
  ENDDO 
  DEALLOCATE(fi,fo,rfmt1)

  RETURN 
END SUBROUTINE 

