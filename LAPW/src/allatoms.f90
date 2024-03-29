SUBROUTINE allatoms()
  USE m_constants, ONLY: solsc
  USE m_atoms, ONLY: nspecies
  USE m_atomic_species, ONLY: xctsp, nstspmax, nstsp, nsp, nrspmax, lsp, vrsp, &
       ksp, rhosp, occsp, rsp, evalsp, ptnucl, spzn, nrsp, spfname
  USE modxcifc, ONLY: getxcdata
  IMPLICIT NONE 
  LOGICAL :: hybrid_
  INTEGER :: xcspin_,xcgrad_
  INTEGER :: is, ist
  REAL(8) :: hybridc_
  CHARACTER(512) :: xcdescr_
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rwf(:,:,:)

  ! allocate global species charge density and potential arrays
  IF(allocated(rhosp)) DEALLOCATE(rhosp)
  ALLOCATE(rhosp(nrspmax,nspecies))
  IF(allocated(vrsp)) DEALLOCATE(vrsp)
  ALLOCATE(vrsp(nrspmax,nspecies))

  ! get the exchange-correlation functional data
  CALL getxcdata(xctsp,xcdescr_,xcspin_,xcgrad_,hybrid_,hybridc_)

  ALLOCATE(rwf(nrspmax,2,nstspmax))
  DO is=1,nspecies

    CALL atom(solsc,ptnucl,spzn(is),nstsp(is),nsp(:,is),lsp(:,is),ksp(:,is), &
     occsp(:,is),xctsp,xcgrad_,nrsp(is),rsp(:,is),evalsp(:,is),rhosp(:,is), &
     vrsp(:,is),rwf)

    DO ist=1,nstsp(is)
      write(*,'(1x,I3,F18.10)') ist, evalsp(ist,is)
    ENDDO

    write(*,'(1x,A,ES18.10)') 'sum(rhosp) = ', sum(rhosp)

  ENDDO 
  DEALLOCATE(rwf)
  
  RETURN 

END SUBROUTINE 
