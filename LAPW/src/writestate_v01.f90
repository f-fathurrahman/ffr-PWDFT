SUBROUTINE writestate()
use modmain
! !DESCRIPTION:
!   Writes the charge density, potentials and other relevant variables to the
!   file {\tt STATE.OUT}. Note to developers: changes to the way the variables
!   are written should be mirrored in {\tt readstate}.
IMPLICIT NONE 
! local variables
INTEGER :: idm,is,ias
! ALLOCATABLE arrays
REAL(8), ALLOCATABLE :: rfmt(:,:,:),rvfmt(:,:,:,:),rvfcmt(:,:,:,:)
open(40,file='STATE'//trim(filext),form='UNFORMATTED')
WRITE(40) version
WRITE(40) spinpol
WRITE(40) nspecies
WRITE(40) lmmaxo
WRITE(40) nrmtmax
WRITE(40) nrcmtmax
DO is=1,nspecies
  WRITE(40) natoms(is)
  WRITE(40) nrmt(is)
  WRITE(40) rsp(1:nrmt(is),is)
  WRITE(40) nrcmt(is)
  WRITE(40) rcmt(1:nrcmt(is),is)
ENDDO 
WRITE(40) ngridg
WRITE(40) ngvec
WRITE(40) ndmag
WRITE(40) nspinor
WRITE(40) fsmtype
WRITE(40) xcgrad
! muffin-tin functions are unpacked to maintain backward compatibility
ALLOCATE(rfmt(lmmaxo,nrmtmax,natmtot))
IF(spinpol) THEN 
  ALLOCATE(rvfmt(lmmaxo,nrmtmax,natmtot,ndmag))
  ALLOCATE(rvfcmt(lmmaxo,nrcmtmax,natmtot,ndmag))
ENDIF 
! write the density
DO ias=1,natmtot
  is=idxis(ias)
  CALL rfmtpack(.false.,nrmt(is),nrmti(is),rhomt(:,ias),rfmt(:,:,ias))
ENDDO 
WRITE(40) rfmt,rhoir
! write the Coulomb potential
DO ias=1,natmtot
  is=idxis(ias)
  CALL rfmtpack(.false.,nrmt(is),nrmti(is),vclmt(:,ias),rfmt(:,:,ias))
ENDDO 
WRITE(40) rfmt,vclir
! write the exchange-correlation potential
DO ias=1,natmtot
  is=idxis(ias)
  CALL rfmtpack(.false.,nrmt(is),nrmti(is),vxcmt(:,ias),rfmt(:,:,ias))
ENDDO 
WRITE(40) rfmt,vxcir
! write the Kohn-Sham effective potential
DO ias=1,natmtot
  is=idxis(ias)
  CALL rfmtpack(.false.,nrmt(is),nrmti(is),vsmt(:,ias),rfmt(:,:,ias))
ENDDO 
WRITE(40) rfmt,vsir
IF(spinpol) THEN 
! write the magnetisation, exchange-correlation and effective magnetic fields
  DO idm=1,ndmag
    DO ias=1,natmtot
      is=idxis(ias)
      CALL rfmtpack(.false.,nrmt(is),nrmti(is),magmt(:,ias,idm), &
       rvfmt(:,:,ias,idm))
    ENDDO 
  ENDDO 
  WRITE(40) rvfmt,magir
  DO idm=1,ndmag
    DO ias=1,natmtot
      is=idxis(ias)
      CALL rfmtpack(.false.,nrmt(is),nrmti(is),bxcmt(:,ias,idm), &
       rvfmt(:,:,ias,idm))
    ENDDO 
  ENDDO 
  WRITE(40) rvfmt,bxcir
  DO idm=1,ndmag
    DO ias=1,natmtot
      is=idxis(ias)
      CALL rfmtpack(.false.,nrcmt(is),nrcmti(is),bsmt(:,ias,idm), &
       rvfcmt(:,:,ias,idm))
    ENDDO 
  ENDDO 
  WRITE(40) rvfcmt,bsir
! write fixed spin moment magnetic fields
  IF(fsmtype.ne.0) THEN 
    WRITE(40) bfsmc
    WRITE(40) bfsmcmt
  ENDIF 
ENDIF 
! write the tau-DFT exchange-correlation potential
IF(xcgrad.eq.4) THEN 
  DO ias=1,natmtot
    is=idxis(ias)
    CALL rfmtpack(.false.,nrmt(is),nrmti(is),wxcmt(:,ias),rfmt(:,:,ias))
  ENDDO 
  WRITE(40) rfmt,wxcir
ENDIF 
close(40)
DEALLOCATE(rfmt)
IF(spinpol) DEALLOCATE(rvfmt,rvfcmt)
RETURN 
END SUBROUTINE 
!EOC

