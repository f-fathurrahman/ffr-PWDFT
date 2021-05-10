SUBROUTINE allatoms()

use modmain, only: nspecies, xctsp, nstspmax, nstsp, nsp, nrspmax, lsp, vrsp, &
                   ksp, rhosp, occsp, rsp, evalsp, solsc, ptnucl, spzn, nrsp, spfname
use modxcifc, only: getxcdata

  IMPLICIT NONE 
  logical :: hybrid_
  INTEGER :: xcspin_,xcgrad_
  INTEGER :: is
  REAL(8) :: hybridc_
  character(512) :: xcdescr_
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
  ENDDO 
  DEALLOCATE(rwf)
  
  RETURN 

END SUBROUTINE 
