module m_electric_vector_pot

!---------------------------------------------!
!     electric field and vector potential     !
!---------------------------------------------!
! tefield is .true. if a polarising constant electric field is applied
logical tefield
! electric field vector in Cartesian coordinates
real(8) efieldc(3)
! electric field vector in lattice coordinates
real(8) efieldl(3)
! tafield is .true. if a constant vector potential is applied
logical tafield
! vector potential A-field which couples to paramagnetic current
real(8) afieldc(3)
! A-field in lattice coordinates
real(8) afieldl(3)

end module