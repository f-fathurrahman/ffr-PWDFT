module m_lattice

!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8) avec(3,3)
! magnitude of random displacements added to lattice vectors
real(8) rndavec
! inverse of lattice vector matrix
real(8) ainv(3,3)
! reciprocal lattice vectors
real(8) bvec(3,3)
! inverse of reciprocal lattice vector matrix
real(8) binv(3,3)
! unit cell volume
real(8) omega
! Brillouin zone volume
real(8) omegabz
! any vector with length less than epslat is considered zero
real(8) epslat

end module
