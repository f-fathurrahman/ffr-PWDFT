!----------------------------------
SUBROUTINE init_crystal_structure()
!----------------------------------
  USE m_constants, ONLY: fourpi
  USE m_atoms, ONLY: atposl, natoms, nspecies, atposc
  USE m_lattice, ONLY: avec, ainv, binv, bvec, omegabz, omega, epslat
  USE m_muffin_tins, ONLY: rmt, omegamt
  USE m_spin, ONLY: vqlss, vqcss
  USE m_dos_optics_response, ONLY: vecql, vecqc
  IMPLICIT NONE 
  INTEGER :: ia, is

  WRITE(*,*) 'Setting up crystal structure'

  !----------------------------------!
  !     crystal structure set up     !
  !----------------------------------!
  
  ! generate the reciprocal lattice vectors and unit cell volume
  CALL reciplat(avec, bvec, omega, omegabz)
  
  ! inverse of the lattice vector matrix
  CALL r3minv(avec, ainv)
  
  ! inverse of the reciprocal vector matrix
  CALL r3minv(bvec, binv)
  
  ! Cartesian coordinates of the spin-spiral vector
  CALL r3mv(bvec, vqlss, vqcss)
  
  DO is = 1,nspecies
    DO ia = 1,natoms(is)
      ! map atomic lattice coordinates to [0,1)
      CALL r3frac( epslat, atposl(:,ia,is) )
      ! determine atomic Cartesian coordinates
      CALL r3mv( avec, atposl(:,ia,is), atposc(:,ia,is))
    ENDDO 
  ENDDO 
  
  ! check for overlapping muffin-tins and adjust radii if required
  CALL checkmt()
  
  ! compute the total muffin-tin volume (M. Meinert)
  omegamt = 0.d0
  DO is = 1,nspecies
    omegamt = omegamt + dble(natoms(is))*(fourpi/3.d0)*rmt(is)**3
  ENDDO 
  
  ! input q-vector in Cartesian coordinates
  CALL r3mv(bvec, vecql, vecqc)

END SUBROUTINE 
