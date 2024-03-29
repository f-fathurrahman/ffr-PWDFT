module m_mixing

!--------------------------!
!     mixing variables     !
!--------------------------!
! type of mixing to use for the potential
integer mixtype
! mixing type description
character(256) mixdescr
! adaptive mixing parameters (formerly beta0 and betamax)
real(8) amixpm(2)
! subspace dimension for Broyden mixing
integer mixsdb
! Broyden mixing parameters alpha and w0
real(8) broydpm(2)

end module 
