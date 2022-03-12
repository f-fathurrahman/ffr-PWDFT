subroutine write_rhosp()
  USE m_atoms, ONLY: nspecies
  USE m_atomic_species, ONLY: rhosp, rsp, nrsp, spname
  implicit none
  integer :: is, ir
  character(256) :: filename

  write(*,*)
  write(*,*) 'Enter write_rhosp'
  DO is=1,nspecies
    filename = 'rhosp_'//trim(spname(is))//'.dat'
    open(unit=201, file=trim(filename))
    DO ir=1,nrsp(is)
      write(201,*) rsp(ir,is), rhosp(ir,is)
    ENDDO 
    close(201)
    write(*,*) 'rhosp written to file: ', trim(filename)
  ENDDO 
  write(*,*) 'shape of rhosp = ', shape(rhosp)
  write(*,*) 'Exit write_rhosp'
  write(*,*)
end subroutine
