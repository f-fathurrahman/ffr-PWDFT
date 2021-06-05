SUBROUTINE rfplot(np,vpl,rfmt,rfir,fp)
  USE m_atoms, ONLY: natoms, natmtot, idxas, idxis, rsp, nspecies, atposc
  USE m_gvectors, ONLY: ngtot, vgc, igfft, ngvec, ngridg
  USE m_lattice, ONLY: epslat, avec
  USE m_muffin_tins, ONLY: npmtmax, lmmaxo, nrmt, nrmti, nrmtmax, rmt, lmaxo, lmaxi
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: np
  REAL(8), intent(in) :: vpl(3,np)
  REAL(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
  REAL(8), intent(out) :: fp(np)
  ! local variables
  INTEGER ias,is,ip
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt1(:,:,:)
  COMPLEX(8), ALLOCATABLE :: zfft(:)
  
  ! unpack the muffin-tin function
  ALLOCATE(rfmt1(lmmaxo,nrmtmax,natmtot))
  DO ias=1,natmtot
    is=idxis(ias)
    CALL rfmtpack(.false.,nrmt(is),nrmti(is),rfmt(:,ias),rfmt1(:,:,ias))
  ENDDO 

  ! Fourier transform rfir to G-space
  ALLOCATE(zfft(ngtot))
  zfft(:)=rfir(:)
  CALL zfftifc(3,ngridg,-1,zfft)

  ! begin loop over all points
  DO ip=1,np
    CALL rfip(ip)
  ENDDO 
  DEALLOCATE(rfmt1,zfft)
  RETURN 

END SUBROUTINE 

