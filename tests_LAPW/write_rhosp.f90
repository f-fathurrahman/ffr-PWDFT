subroutine write_rhosp()
  USE m_atoms, ONLY: nspecies
  USE m_atomic_species, ONLY: rhosp, rsp, nrsp, spname
  implicit none
  integer :: is, ir

  DO is=1,nspecies
    open(unit=201, file='rhosp_'//trim(spname(is))//'.dat')
    DO ir=1,nrsp(is)
      write(201,*) rsp(ir,is), rhosp(ir,is)
    ENDDO 
    close(201)
  ENDDO 
  write(*,*) 'shape of rhosp = ', shape(rhosp)
end subroutine
