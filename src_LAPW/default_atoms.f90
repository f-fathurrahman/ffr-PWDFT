SUBROUTINE default_atoms()
  USE m_atoms, ONLY: natoms, xctsp, atposl, atposc, ecvcut, esccut, molecule, nspecies, &
               primcell, ptnucl, rndatposc, sppath
  IMPLICIT NONE 

  natoms(:) = 0
  xctsp(1) = 3
  xctsp(2:3) = 0
  atposl(:,:,:) = 0.d0
  atposc(:,:,:) = 0.d0
  ecvcut = -3.5d0
  molecule=.false.
  nspecies=0
  primcell=.false.
  ptnucl=.true.
  rndatposc=0.d0
  sppath=''
  esccut=-0.4d0

END SUBROUTINE 