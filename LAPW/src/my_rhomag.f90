!---------------------
SUBROUTINE my_rhomag()
!---------------------
  USE m_atoms, ONLY: natmtot, idxis
  USE m_spin, ONLY: spinpol, nspnfv, ndmag, ncmag
  USE m_muffin_tins, ONLY: npcmt, nrcmti, nrcmt, npmtmax, lmmaxapw
  USE m_density_pot_xc, ONLY: rhomt, rhoir, magir, magmt
  USE m_apwlo, ONLY: apwordmax
  USE m_states, ONLY: nstfv, nstsv, occsv
  USE m_kpoints, ONLY: nkpt, vkl, wkpt
  USE m_gvectors, ONLY: ngtot
  USE m_gkvectors, ONLY: ngkmax, ngk, igkig, vgkc, vgkl, gkc, sfacgk
  USE m_hamiltonian, ONLY: nmatmax
  USE m_misc, ONLY: filext
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,ispn,idm
  INTEGER :: is,ias,n
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)

  write(*,*)
  write(*,*) '---------------'
  write(*,*) 'ENTER my_rhomag'
  write(*,*) '---------------'
  write(*,*)


  ! set the charge density and magnetisation to zero
  !
  ! rho, muffin tin
  DO ias = 1,natmtot
    is = idxis(ias)
    rhomt(1:npcmt(is),ias) = 0.d0
  ENDDO 
  !
  ! interstitial
  rhoir(:) = 0.d0
  !
  ! magnetization, muffin tin
  DO idm = 1,ndmag
    DO ias = 1,natmtot
      is = idxis(ias)
      magmt(1:npcmt(is),ias,idm)=0.d0
    ENDDO
    ! interstitial
    magir(:,idm) = 0.d0
  ENDDO 

  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))

  DO ik=1,nkpt
    ! get the eigenvectors from file
    CALL my_getevecfv(filext, ik, vkl(:,ik), vgkl(:,:,:,ik),evecfv)
    CALL getevecsv(filext, ik, vkl(:,ik), evecsv)
    ! find the matching coefficients
    DO ispn=1,nspnfv
      CALL match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    ENDDO 
    ! add to the density and magnetisation
    ! pass apwalm and eigenvectors
    CALL my_rhomagk(ngk(:,ik),igkig(:,:,ik),wkpt(ik),occsv(:,ik),apwalm, &
     evecfv,evecsv)
  ENDDO 

  DEALLOCATE(apwalm,evecfv,evecsv)

  ! convert muffin-tin density/magnetisation to spherical harmonics
  CALL rhomagsh()

  ! symmetrise the density
  CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,rhomt,rhoir)

  ! convert the density from a coarse to a fine radial mesh
  CALL rfmtctof(rhomt)

  ! symmetrise the magnetisation
  IF(spinpol) CALL symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,magmt,magir)

  ! convert the magnetisation from a coarse to a fine radial mesh
  DO idm=1,ndmag
    CALL rfmtctof(magmt(:,:,idm))
  ENDDO 

  ! add the core density and magnetisation to the total
  CALL rhocore()

  ! calculate the charges
  CALL charge()

  ! normalise the density
  CALL rhonorm()

  ! calculate the moments
  IF(spinpol) CALL moment()

  write(*,*)
  write(*,*) '---------------'
  write(*,*) 'EXIT my_rhomag'
  write(*,*) '---------------'
  write(*,*)

  RETURN 
END SUBROUTINE 
