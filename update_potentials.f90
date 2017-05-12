SUBROUTINE update_potentials()

  USE m_realspace, ONLY : Npoints
  USE m_hamiltonian, ONLY : Rhoe, V_Hartree, V_xc

  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  
  ALLOCATE( epsxc(Npoints) )
  ALLOCATE( depsxc(Npoints) )

  CALL Poisson_solve_fft( Rhoe, V_Hartree )

  CALL excVWN( Npoints, Rhoe, epsxc )
  CALL excpVWN( Npoints, Rhoe, depsxc )

  V_xc(:) = epsxc(:) + Rhoe(:)*depsxc(:)

  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END SUBROUTINE 
