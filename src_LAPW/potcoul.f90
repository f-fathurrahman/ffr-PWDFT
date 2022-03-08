SUBROUTINE potcoul()
  USE m_atoms, ONLY: idxis, natmtot
  USE m_atomic_species, ONLY: vcln
  USE m_muffin_tins, ONLY: nrmt, nrmti, nrmtmax, npmtmax, npmti, wprmt, rlmt, npmt, &
            lmmaxi, lmmaxo
  USE m_gvectors, ONLY: ylmg, sfacg, ngtot, ngvec, ngridg, jlgrmt, igfft, gclg, gc
  USE m_electric_vector_pot, ONLY: tefield
  USE m_density_pot_xc, ONLY: vclir, rhoir, rhomt, vclmt
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ias
  INTEGER :: nr,nri,ir,i
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: zrhomt(:,:),zrhoir(:)
  COMPLEX(8), ALLOCATABLE :: zvclmt(:,:),zvclir(:)
  
  ALLOCATE(zrhomt(npmtmax,natmtot))
  
  ! convert real muffin-tin charge density to complex spherical harmonic expansion
  DO ias=1,natmtot
    is=idxis(ias)
    CALL rtozfmt(nrmt(is),nrmti(is),rhomt(:,ias),zrhomt(:,ias))
  ENDDO 
  
  ! solve the complex Poisson's equation in the muffin-tins
  ALLOCATE(zvclmt(npmtmax,natmtot))
  CALL genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,zrhomt,zvclmt)
  DEALLOCATE(zrhomt)
  
  ! add the nuclear monopole potentials
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    i=1
    DO ir=1,nri
      zvclmt(i,ias) = zvclmt(i,ias) + vcln(ir,is)
      i = i + lmmaxi
    ENDDO 
    DO ir=nri+1,nr
      zvclmt(i,ias) = zvclmt(i,ias)+vcln(ir,is)
      i=i+lmmaxo
    ENDDO 
  ENDDO 
  
  ! store real interstitial charge density in complex array
  ALLOCATE(zrhoir(ngtot))
  zrhoir(:)=rhoir(:)
  ! solve Poisson's equation in the entire unit cell
  ALLOCATE(zvclir(ngtot))
  CALL zpotcoul(nrmt,nrmti,npmt,npmti,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
   ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)

  ! convert complex muffin-tin potential to real spherical harmonic expansion
  DO ias=1,natmtot
    is=idxis(ias)
    CALL ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),vclmt(:,ias))
  ENDDO 
  
  ! store complex interstitial potential in real array
  vclir(:)=dble(zvclir(:))
  DEALLOCATE(zrhoir,zvclmt,zvclir)

  ! apply constant electric field if required
  IF(tefield) CALL potefield()
  
  RETURN 

END SUBROUTINE 