! !INPUT/OUTPUT PARAMETERS:
!   sol    : speed of light in atomic units (in,real)
!   ptnucl : .true. if the nucleus is a point particle (in,logical)
!   zn     : nuclear charge (in,real)
!   nst    : number of states to solve for (in,integer)
!   n      : priciple quantum number of each state (in,integer(nst))
!   l      : quantum number l of each state (in,integer(nst))
!   k      : quantum number k (l or l+1) of each state (in,integer(nst))
!   occ    : occupancy of each state (inout,real(nst))
!   xctype : exchange-correlation type (in,integer(3))
!   xcgrad : 1 for GGA functional, 0 otherwise (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   eval   : eigenvalue without rest-mass energy for each state (out,real(nst))
!   rho    : charge density (out,real(nr))
!   vr     : self-constistent potential (out,real(nr))
!   rwf    : major and minor components of radial wavefunctions for each state
!            (out,real(nr,2,nst))
! !DESCRIPTION:
!   Solves the Dirac-Kohn-Sham equations for an atom using the
!   exchange-correlation functional {\tt xctype} and RETURN s the self-consistent
!   radial wavefunctions, eigenvalues, charge densities and potentials. Requires
!   the exchange-correlation interface routine {\tt xcifc}.
SUBROUTINE atom(sol,ptnucl,zn,nst,n,l,k,occ,xctype,xcgrad,nr,r,eval,rho,vr,rwf)

  USE modxcifc, only: xcifc

  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: sol
  LOGICAL, intent(in) :: ptnucl
  REAL(8), intent(in) :: zn
  INTEGER, intent(in) :: nst
  INTEGER, intent(in) :: n(nst),l(nst),k(nst)
  REAL(8), intent(inout) :: occ(nst)
  INTEGER, intent(in) :: xctype(3),xcgrad
  INTEGER, intent(in) :: nr
  REAL(8), intent(in) :: r(nr)
  REAL(8), intent(out) :: eval(nst)
  REAL(8), intent(out) :: rho(nr),vr(nr)
  REAL(8), intent(out) :: rwf(nr,2,nst)
  ! local variables
  INTEGER, parameter :: maxscl=200
  INTEGER ir,ist,iscl
  REAL(8), parameter :: fourpi=12.566370614359172954d0
  ! potential convergence tolerance
  REAL(8), parameter :: eps=1.d-6
  REAL(8) sum,dv,dvp,ze,beta,t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: vn(:),vh(:),ex(:),ec(:),vx(:),vc(:),vrp(:)
  REAL(8), ALLOCATABLE :: ri(:),wpr(:,:),fr1(:),fr2(:),gr1(:),gr2(:)
  REAL(8), ALLOCATABLE :: grho(:),g2rho(:),g3rho(:)

  IF(nst <= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(atom): invalid nst : ",I8)') nst
    WRITE(*,*)
    stop
  ENDIF 

  WRITE(*,*) '--- from atom ---'
  WRITE(*,*) 'sol    = ', sol
  WRITE(*,*) 'ptnucl = ', ptnucl
  WRITE(*,*) 'zn     = ', zn
  WRITE(*,*) 'nst    = ', nst
  DO ist = 1,nst
    WRITE(*,'(1x,A,I5)', advance="no") ' n = ', n(ist)
    WRITE(*,'(1x,A,I5)', advance="no") ' l = ', l(ist)
    WRITE(*,'(1x,A,I5)', advance="no") ' k = ', k(ist)
    WRITE(*,*)
  enddo
  !k,occ,

  ! allocate local arrays
  ALLOCATE(vn(nr),vh(nr),ex(nr),ec(nr),vx(nr),vc(nr),vrp(nr))
  ALLOCATE(ri(nr),wpr(4,nr),fr1(nr),fr2(nr),gr1(nr),gr2(nr))
  IF(xcgrad==1) THEN 
    ALLOCATE( grho(nr), g2rho(nr), g3rho(nr) )
  ENDIF 

  ! find total electronic charge
  ze = 0.d0
  DO ist = 1,nst
    ze = ze + occ(ist)
  ENDDO 

  ! set up nuclear potential
  WRITE(*,*) 'ptnucl = ', ptnucl
  CALL potnucl(ptnucl,nr,r,zn,vn)

  DO ir=1,nr
    ri(ir)=1.d0/r(ir)
    ! initialise the Kohn-Sham potential to the nuclear potential
    vr(ir)=vn(ir)
  ENDDO 

  ! determine the weights for radial integration
  CALL wsplintp(nr,r,wpr)
  dvp=0.d0
  vrp(:)=0.d0

  ! initialise mixing parameter
  beta=0.5d0
  
  ! initialise eigenvalues to relativistic values (minus the rest mass energy)
  DO ist=1,nst
    t1=sqrt(dble(k(ist)**2)-(zn/sol)**2)
    t1=(dble(n(ist)-abs(k(ist)))+t1)**2
    t1=1.d0+((zn/sol)**2)/t1
    eval(ist)=sol**2/sqrt(t1)-sol**2
  ENDDO 

  ! start of self-consistent loop
  DO iscl=1,maxscl
  
    ! solve the Dirac equation for each state
    DO ist=1,nst
      CALL rdirac(sol, n(ist), l(ist), k(ist), nr, r, vr, &
                  eval(ist), rwf(:,1,ist), rwf(:,2,ist))
    ENDDO 
  
    ! compute the charge density
    DO ir=1,nr
      sum=0.d0
      DO ist=1,nst
        sum=sum+occ(ist)*(rwf(ir,1,ist)**2+rwf(ir,2,ist)**2)
      ENDDO 
      fr1(ir)=sum
      fr2(ir)=sum*ri(ir)
      rho(ir)=(1.d0/fourpi)*sum*ri(ir)**2
    ENDDO 
    CALL splintwp(nr,wpr,fr1,gr1)
    CALL splintwp(nr,wpr,fr2,gr2)
    
    ! find the Hartree potential
    t1 = gr2(nr)
    DO ir = 1,nr
      vh(ir) = gr1(ir)*ri(ir) + t1 - gr2(ir)
    ENDDO 

    ! normalise charge density and potential
    t1=ze/gr1(nr)
    rho(:)=t1*rho(:)
    vh(:)=t1*vh(:)
    ! compute the exchange-correlation energy and potential
    IF(xcgrad==1) THEN 
      ! GGA functional
      ! |grad rho|
      CALL fderiv(1,nr,r,rho,grho)
      ! grad^2 rho
      CALL fderiv(2,nr,r,rho,g2rho)
      DO ir=1,nr
        g2rho(ir) = g2rho(ir) + 2.d0*ri(ir)*grho(ir)
      ENDDO 
      ! approximate (grad rho).(grad |grad rho|)
      DO ir=1,nr
        g3rho(ir) = grho(ir)*g2rho(ir)
      ENDDO 
      CALL xcifc(xctype,n=nr,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
       ec=ec,vx=vx,vc=vc)
    ELSE 
      ! LDA functional
      CALL xcifc(xctype,n=nr,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
    ENDIF 
    
    ! self-consistent potential
    vr(:)=vh(:)+vx(:)+vc(:)
    
    ! determine change in potential
    sum=0.d0
    DO ir=1,nr
      sum=sum+(vr(ir)-vrp(ir))**2
    ENDDO 
    dv=sqrt(sum)/dble(nr)
    
    IF(iscl.gt.2) THEN 
      ! reduce beta if change in potential is diverging
      IF(dv.gt.dvp) beta=beta*0.8d0
      beta=max(beta,0.01d0)
    ENDIF 
    
    dvp=dv
    DO ir=1,nr
      ! mix old and new potentials
      vr(ir)=(1.d0-beta)*vrp(ir)+beta*vr(ir)
      vrp(ir)=vr(ir)
      ! add nuclear potential
      vr(ir)=vr(ir)+vn(ir)
    ENDDO 
    
    ! check for convergence
    IF((iscl.gt.2).and.(dv.lt.eps)) then
      write(*,*) 'Converged at iscl = ', iscl
      write(*,'(1x,A,ES18.10)') 'dv = ', dv
      goto 10
    endif
  
  ENDDO ! end self-consistent loop
  WRITE(*,*)
  WRITE(*,'("Warning(atom): maximum iterations exceeded")')
  
  10 CONTINUE 
  DEALLOCATE(vn,vh,ex,ec,vx,vc,vrp)
  DEALLOCATE(ri,wpr,fr1,fr2,gr1,gr2)
  
  IF(xcgrad==1) DEALLOCATE(grho,g2rho,g3rho)
  
  RETURN 
END SUBROUTINE 