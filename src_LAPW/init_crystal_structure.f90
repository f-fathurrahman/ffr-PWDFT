!----------------------------------
SUBROUTINE init_crystal_structure()
!----------------------------------
  USE modmain, ONLY: &
               avec, nspecies, atposl, ainv, binv, bvec, fourpi, rmt, natoms, &
               omegamt, omegabz, omega, epslat, atposc, &
               vqlss, vecql, vecqc, vqcss
  IMPLICIT NONE 
  INTEGER :: ia, is

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
      call r3frac( epslat, atposl(:,ia,is) )
      ! determine atomic Cartesian coordinates
      call r3mv( avec, atposl(:,ia,is), atposc(:,ia,is))
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
  
  ! write to VARIABLES.OUT
  !call writevars('avec',nv=9,rva=avec)
  !call writevars('bvec',nv=9,rva=bvec)
  !call writevars('omega',rv=omega)
  !do is=1,nspecies
  !  call writevars('atposl',l=is,nv=3*natoms(is),rva=atposl(:,:,is))
  !end DO

END SUBROUTINE 
