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
  write(*,*) '<div> ENTER my_rhomag'
  write(*,*)


  ! set the charge density and magnetisation to zero
  !
  ! rho, muffin tin
  DO ias = 1,natmtot
    is = idxis(ias)
    !rhomt(1:npcmt(is),ias) = 0.d0 ! XXX only 1:npcmt ? only coarse mesh?
    rhomt(:,ias) = 0.d0 ! XXX debug
  ENDDO
  !
  ! interstitial
  rhoir(:) = 0.d0
  !
  ! magnetization, muffin tin
  DO idm = 1,ndmag
    DO ias = 1,natmtot
      is = idxis(ias)
      !magmt(1:npcmt(is),ias,idm) = 0.d0
      magmt(:,ias,idm) = 0.d0
    ENDDO
    ! interstitial
    magir(:,idm) = 0.d0
  ENDDO

  ALLOCATE( apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv) )
  ALLOCATE( evecfv(nmatmax,nstfv,nspnfv), evecsv(nstsv,nstsv) )

  ! eigenvectors (wavefunctions are used and read in this loop only)
  DO ik=1,nkpt
    ! get the eigenvectors from file
    CALL getevecfv(filext, ik, vkl(:,ik), vgkl(:,:,:,ik),evecfv)
    CALL getevecsv(filext, ik, vkl(:,ik), evecsv)
    ! find the matching coefficients
    ! ffr: do they depend on evecfv ? apwfr ? try debug this first in ElkDFTWrapper
    DO ispn = 1,nspnfv
      CALL match( ngk(ispn,ik), vgkc(:,:,ispn,ik), gkc(:,ispn,ik), &
                  sfacgk(:,:,ispn,ik), apwalm(:,:,:,:,ispn))
    ENDDO
    ! add to the density and magnetisation
    ! pass apwalm and eigenvectors
    CALL my_rhomagk( ik, ngk(:,ik), igkig(:,:,ik), wkpt(ik), occsv(:,ik), apwalm, &
                     evecfv, evecsv )
  ENDDO
  ! matching coefs, eigenvectors are no longer needed
  DEALLOCATE(apwalm, evecfv, evecsv)

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag'
  !RETURN ! DEBUG

  ! convert muffin-tin density/magnetisation to spherical harmonics
  CALL rhomagsh()

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after rhomagsh'
  !RETURN ! DEBUG

  ! symmetrise the density
  write(*,*) 'before symrf, sum(rhomt) = ', sum(rhomt)
  write(*,*) '              sum(rhoir) = ', sum(rhoir)
  CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,rhomt,rhoir)
  write(*,*) 'after symrf, sum(rhomt) = ', sum(rhomt)
  write(*,*) '             sum(rhoir) = ', sum(rhoir)

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after symrf'
  !RETURN ! DEBUG

  ! convert the density from a coarse to a fine radial mesh
  CALL rfmtctof(rhomt)

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after rfmtctof'
  !RETURN ! DEBUG

  ! symmetrise the magnetisation
  IF(spinpol) CALL my_symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,magmt,magir)

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after symrvf'
  !RETURN ! DEBUG

  ! convert the magnetisation from a coarse to a fine radial mesh
  DO idm=1,ndmag
    CALL rfmtctof(magmt(:,:,idm))
  ENDDO

  ! add the core density and magnetisation to the total
  CALL rhocore()

  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after rhocore'
  !RETURN ! DEBUG

  ! calculate the charges
  CALL charge()
  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after charge'
  !RETURN ! DEBUG

  ! normalise the density
  CALL rhonorm()
  !write(*,*)
  !write(*,*) '>>>>> EARLY RETURN in my_rhomag after rhonorm'
  !RETURN ! DEBUG

  ! calculate the moments
  IF(spinpol) CALL moment()

  write(*,*)
  write(*,*) '</div> EXIT my_rhomag'
  write(*,*)

  RETURN
END SUBROUTINE
