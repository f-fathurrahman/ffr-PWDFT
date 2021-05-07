SUBROUTINE genapwlofr()
  IMPLICIT NONE 
  CALL genapwfr() ! generate the APW radial functions
  CALL genlofr() ! generate the local-orbital radial functions
  CALL olprad() ! compute the overlap radial integrals
  CALL hmlrad() ! compute the Hamiltonian radial integrals
  RETURN 
END SUBROUTINE 