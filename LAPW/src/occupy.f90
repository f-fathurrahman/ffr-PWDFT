SUBROUTINE occupy()
  use m_convergence, only: iscl
  ! iscl is used when autoswidth is active, which in this case, only
  ! when iscl > 1
  use m_kpoints, only: wkpt, nkpt
  use m_states, only: bandgap, evalsv, ikgap, occsv, swidth, stype, nstsv, occmax,&
                    & fermidos, efermi, epsocc, autoswidth
  use m_apwlo, only: e0min
  use m_charge_moment_current, only: chgval
  use m_apwlo, only: e0min  
  !
  ! actual inputs: evalsv
  ! actual output: occsv
  !
  ! !DESCRIPTION:
  !   Finds the Fermi energy and sets the occupation numbers for the
  !   second-variational states using the routine {\tt fermi}.
  !
  IMPLICIT NONE 
  ! local variables
  INTEGER, parameter :: maxit=1000
  INTEGER ik,ist,it
  REAL(8) e0,e1,e
  REAL(8) chg,x,t1
  ! external functions
  REAL(8) sdelta,stheta
  external sdelta,stheta

  ! determine the smearing width automatically if required
  IF((autoswidth).and.(iscl.gt.1)) CALL findswidth
  ! find minimum and maximum eigenvalues
  e0=evalsv(1,1)
  e1=e0
  DO ik=1,nkpt
    DO ist=1,nstsv
      e=evalsv(ist,ik)
      IF(e.lt.e0) e0=e
      IF(e.gt.e1) e1=e
    ENDDO 
  ENDDO 
  IF(e0.lt.e0min) THEN 
    WRITE(*,*)
    WRITE(*,'("Warning(occupy): minimum eigenvalue less than minimum &
     &linearisation energy : ",2G18.10)') e0,e0min
    WRITE(*,'(" for s.c. loop ",I5)') iscl
  ENDIF 
  t1=1.d0/swidth
  ! determine the Fermi energy using the bisection method
  DO it=1,maxit
    efermi=0.5d0*(e0+e1)
    chg=0.d0
    DO ik=1,nkpt
      DO ist=1,nstsv
        e=evalsv(ist,ik)
        IF(e.lt.e0min) THEN 
          occsv(ist,ik)=0.d0
        else
          x=(efermi-e)*t1
          occsv(ist,ik)=occmax*stheta(stype,x)
          chg=chg+wkpt(ik)*occsv(ist,ik)
        ENDIF 
      ENDDO 
    ENDDO 
    IF(chg.lt.chgval) THEN 
      e0=efermi
    else
      e1=efermi
    ENDIF 
    IF((e1-e0).lt.1.d-12) goto 10
  ENDDO 
  WRITE(*,*)
  WRITE(*,'("Warning(occupy): could not find Fermi energy")')
  10 continue
  ! find the density of states at the Fermi surface in units of
  ! states/Hartree/unit cell
  fermidos=0.d0
  DO ik=1,nkpt
    DO ist=1,nstsv
      x=(evalsv(ist,ik)-efermi)*t1
      fermidos=fermidos+wkpt(ik)*sdelta(stype,x)*t1
    ENDDO 
    IF(abs(occsv(nstsv,ik)).gt.epsocc) THEN 
      WRITE(*,*)
      WRITE(*,'("Warning(occupy): not enough empty states for k-point ",I6)') ik
      WRITE(*,'(" and s.c. loop ",I5)') iscl
    ENDIF 
  ENDDO 
  fermidos=fermidos*occmax

  ! estimate the indirect band gap (FC)
  e0=-1.d8
  e1=1.d8
  ikgap(1)=1
  ikgap(2)=1
  ! these loops are incorrectly ordered to fix a bug in versions 17 and 18 of the
  ! Intel compiler
  DO ist=1,nstsv
    DO ik=1,nkpt
      e=evalsv(ist,ik)
      IF(e.lt.efermi) THEN 
        IF(e.gt.e0) THEN 
          e0=e
          ikgap(1)=ik
        ENDIF 
      else
        IF(e.lt.e1) THEN 
          e1=e
          ikgap(2)=ik
        ENDIF 
      ENDIF 
    ENDDO 
  ENDDO 
  bandgap(1)=e1-e0

  ! write band gap to test file
  WRITE(*,*) 'estimated indirect band gap: ', bandgap(1)

  ! estimate the direct band gap
  e=1.d8
  ikgap(3)=1
  DO ik=1,nkpt
    e0=-1.d8
    e1=1.d8
    DO ist=1,nstsv
      t1=evalsv(ist,ik)
      IF(t1.le.efermi) THEN 
        IF(t1.gt.e0) e0=t1
      else
        IF(t1.lt.e1) e1=t1
      ENDIF 
    ENDDO 
    t1=e1-e0
    IF(t1.lt.e) THEN 
      e=t1
      ikgap(3)=ik
    ENDIF 
  ENDDO 
  bandgap(2)=e
  RETURN 
END SUBROUTINE 
