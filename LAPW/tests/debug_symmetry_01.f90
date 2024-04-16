program debug_main
  IMPLICIT NONE 

  ! read input and initialize some variables
  CALL read_input()
  CALL init_timing()

  CALL init_am_variables()
  CALL init_idx_atom_species()
  CALL init_spin_variables()
  
  CALL init_crystal_structure()
  
  CALL init_vector_field_E_A()
  
  CALL symmetry()

end program

