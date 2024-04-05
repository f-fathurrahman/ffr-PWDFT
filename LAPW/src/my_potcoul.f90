!----------------------
SUBROUTINE my_potcoul()
!----------------------
  USE m_atoms, ONLY: idxis, natmtot
  USE m_atomic_species, ONLY: vcln
  USE m_muffin_tins, ONLY: nrmt, nrmti, nrmtmax, npmtmax, npmti, wprmt, rlmt, npmt, &
            lmmaxi, lmmaxo
  USE m_gvectors, ONLY: ylmg, sfacg, ngtot, ngvec, ngridg, jlgrmt, igfft, gclg, gc
  USE m_density_pot_xc, ONLY: rhoir, rhomt, vclir, vclmt
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ias
  INTEGER :: nr,nri,ir,i
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: zrhomt(:,:), zrhoir(:)
  COMPLEX(8), ALLOCATABLE :: zvclmt(:,:), zvclir(:)

  write(*,*)  
  write(*,*) 'Enter my_potcoul'

  ALLOCATE(zrhomt(npmtmax,natmtot))
  
  ! convert real muffin-tin charge density to complex spherical harmonic expansion
  DO ias=1,natmtot
    is = idxis(ias)
    CALL r_to_zf_mt(nrmt(is), nrmti(is), rhomt(:,ias), zrhomt(:,ias))
  ENDDO 
  write(*,*) 'shape rhomt = ', shape(rhomt)
  write(*,*) 'shape zrhomt = ', shape(zrhomt)

  write(*,*) 'sum rhomt = ', sum(rhomt)
  write(*,*) 'sum zrhomt = ', sum(zrhomt)  

  ! solve the complex Poisson's equation in the muffin-tins
  ALLOCATE(zvclmt(npmtmax,natmtot))
  CALL my_genzvclmt(nrmt, nrmti, nrmtmax, rlmt, wprmt, npmtmax, zrhomt, zvclmt)
  DEALLOCATE(zrhomt)

  write(*,*) 'shape(zvclmt) = ', shape(zvclmt)
  write(*,*) 'sum(zvclmt) = ', sum(zvclmt)

  ! add the nuclear monopole potentials
  DO ias = 1,natmtot
    is = idxis(ias)
    nr = nrmt(is)
    nri = nrmti(is)
    i = 1
    DO ir = 1,nri
      zvclmt(i,ias) = zvclmt(i,ias) + vcln(ir,is)
      i = i + lmmaxi
    ENDDO 
    DO ir = nri+1, nr
      zvclmt(i,ias) = zvclmt(i,ias) + vcln(ir,is)
      i = i + lmmaxo
    ENDDO 
  ENDDO 

  write(*,*) 'sum(vcln) = ', sum(vcln)
  write(*,*) 'After adding vcln: sum(zvclmt) = ', sum(zvclmt)
  
  ! store real interstitial charge density in complex array
  ALLOCATE(zrhoir(ngtot))
  zrhoir(:) = rhoir(:)
  write(*,*) 'In potcoul, sum rhoir = ', sum(rhoir)
  write(*,*) 'In potcoul, sum zrhoir = ', sum(zrhoir)

  ! solve Poisson's equation in the entire unit cell
  ALLOCATE(zvclir(ngtot))
  CALL my_zpotcoul(nrmt,nrmti,npmt,npmti,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
   ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)

  ! convert complex muffin-tin potential to real spherical harmonic expansion
  vclmt(:,:)  = 0.d0 ! ffr
  DO ias=1,natmtot
    is=idxis(ias)
    CALL z_to_rf_mt(nrmt(is),nrmti(is),zvclmt(:,ias),vclmt(:,ias))
  ENDDO 
  
  ! store complex interstitial potential in real array
  vclir(:) = dble(zvclir(:))

  write(*,*) 'sum vclmt = ', sum(vclmt)
  write(*,*) 'sum abs vclir = ', sum(abs(vclir))  

  DEALLOCATE(zrhoir, zvclmt, zvclir)

  write(*,*)  
  write(*,*) 'Exit my_potcoul'

  RETURN 

END SUBROUTINE 