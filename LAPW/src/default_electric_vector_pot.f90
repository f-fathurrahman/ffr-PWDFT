SUBROUTINE default_electric_vector_pot()
  USE m_electric_vector_pot, ONLY: afieldc, efieldc
  efieldc(:)=0.d0
  afieldc(:)=0.d0
END SUBROUTINE