subroutine write_radial_mt()
  use m_muffin_tins
  implicit none

  write(*,*)
  write(*,*) 'nrmt = ', nrmt
  write(*,*) 'shape(rlmt) = ', shape(rlmt)
  write(99,*) rlmt(1:nrmt(1),1,1)
end subroutine