PROGRAM lapwdft
  CALL read_input()

  CALL init_timing()
  CALL init_am_variables()
  CALL init_idx_atom_species()
  CALL init_spin_variables()
  CALL init_crystal_structure()
  CALL init_vector_field_E_A()
  CALL symmetry() ! crystal symmetry set up
  CALL init_radial_meshes()
  CALL init_charges_states()
  CALL init_gvector_arrays()
  CALL init_atoms_cores()
  CALL init_chgden_pot_xc()
  CALL init_forces()

  !CALL info_lattice()
  !CALL info_apwlo()

END PROGRAM