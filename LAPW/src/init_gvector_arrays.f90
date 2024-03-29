!-------------------------------
SUBROUTINE init_gvector_arrays()
!-------------------------------

  USE m_atoms, ONLY: nspecies, natmtot, natoms
  USE m_muffin_tins, ONLY: rmt, lmaxo
  USE m_density_pot_xc, ONLY: npsd, lnpsd
  USE m_gvectors, ONLY: vgc, sfacg, ngvec, gmaxvr, jlgrmt, gc, ffacg, ngtot
  USE m_gkvectors, ONLY: isgkmax, gkmax, rgkmax

  IMPLICIT NONE 
  INTEGER :: is
  REAL(8) :: t1, rsum

  WRITE(*,*) 'Setting up Gvector arrays'

  !-------------------------!
  !     G-vector arrays     !
  !-------------------------!
  ! determine gkmax from rgkmax
  IF (nspecies.eq.0) isgkmax=-2
  select case(isgkmax)
  case(:-4)
    write(*,*) 'Use largest muffin-tin radius'
    gkmax=rgkmax/maxval(rmt(1:nspecies))
  case(-3)
    write(*,*) 'Use smallest muffin-tin radius'
    gkmax=rgkmax/minval(rmt(1:nspecies))
  case(-2)
    write(*,*) 'Use the fixed value of 2.0'
    gkmax=rgkmax/2.d0
  case(-1)
    write(*,*) 'Use average muffin-tin radius'
    rsum=0.d0
    DO is=1,nspecies
      write(*,*) 'rmt(is) = ', rmt(is)
      rsum=rsum+dble(natoms(is))*rmt(is)
    ENDDO 
    rsum=rsum/dble(natmtot)
    gkmax=rgkmax/rsum
    
    write(*,*) 'rgkmax = ', rgkmax
    write(*,*) 'rsum   = ', rsum
    write(*,*) 'gkmax  = ', gkmax

  case(1:)
    write(*,*) 'Use user-specified muffin-tin radius'
    IF(isgkmax <= nspecies) THEN 
      gkmax = rgkmax/rmt(isgkmax)
    ELSE 
      WRITE(*,*)
      WRITE(*,'("Error(init0): isgkmax > nspecies : ",2I8)') isgkmax,nspecies
      WRITE(*,*)
      STOP 
    ENDIF 
  END SELECT 

  ! generate the G-vectors
  CALL gengvec()
  
  ! Poisson solver pseudocharge density constant
  IF(nspecies.gt.0) THEN 
    t1=0.25d0*gmaxvr*maxval(rmt(1:nspecies))
  ELSE 
    t1=0.25d0*gmaxvr*2.d0
  ENDIF 
  npsd=max(nint(t1),1)
  lnpsd=lmaxo+npsd+1
  
  ! generate the Coulomb Green's function in G-space = fourpi / G^2
  CALL gengclg()
  
  ! compute the spherical Bessel functions j_l(|G|R_mt)
  IF(allocated(jlgrmt)) DEALLOCATE(jlgrmt)
  ALLOCATE(jlgrmt(0:lnpsd,ngvec,nspecies))
  CALL genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)
  
  ! generate the spherical harmonics of the G-vectors
  CALL genylmg()
  
  ! allocate structure factor array for G-vectors
  IF(allocated(sfacg)) DEALLOCATE(sfacg)
  ALLOCATE(sfacg(ngvec,natmtot))
  
  ! generate structure factors for G-vectors
  CALL gensfacgp(ngvec,vgc,ngvec,sfacg)
  
  ! generate the smooth step function form factors
  IF(allocated(ffacg)) DEALLOCATE(ffacg)
  ALLOCATE(ffacg(ngtot,nspecies))
  DO is=1,nspecies
    CALL genffacgp(is,gc,ffacg(:,is))
  ENDDO 
  
  ! generate the characteristic function
  CALL gencfun()

END SUBROUTINE 

