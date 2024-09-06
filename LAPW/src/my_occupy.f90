!---------------------
SUBROUTINE my_occupy()
!---------------------
  use m_convergence, only: iscl
  ! iscl is used when autoswidth is active, which in this case, only
  ! when iscl > 1
  use m_kpoints, only: wkpt, nkpt
  use m_states, only: bandgap, evalsv, ikgap, occsv, swidth, stype, nstsv, occmax,&
                    & fermidos, efermi, epsocc, autoswidth
  ! actual inputs: evalsv
  ! actual output: occsv
  !
  use m_charge_moment_current, only: chgval
  use m_apwlo, only: e0min
  IMPLICIT NONE 
  ! local variables
  INTEGER, parameter :: maxit=1000
  INTEGER ik,ist,it
  REAL(8) e0,e1,e
  REAL(8) chg,x,t1
  ! external functions
  REAL(8) sdelta,stheta
  external sdelta,stheta

  ! ffr: why use evalsv? (second variational)



  ! determine the smearing width automatically if required
  IF( (autoswidth) .and. (iscl > 1)) CALL findswidth()
  !
  ! find minimum and maximum eigenvalues
  e0 = evalsv(1,1)
  e1 = e0
  DO ik = 1,nkpt
    DO ist = 1,nstsv
      e = evalsv(ist,ik)
      IF(e < e0) e0 = e
      IF(e > e1) e1 = e
    ENDDO 
  ENDDO 
  IF(e0 < e0min) THEN 
    WRITE(*,*)
    WRITE(*,'("Warning(my_occupy): minimum eigenvalue less than minimum &
     &linearisation energy : ",2G18.10)') e0,e0min
    WRITE(*,'(" for s.c. loop ",I5)') iscl
  ENDIF 
  t1 = 1.d0/swidth
  ! determine the Fermi energy using the bisection method
  DO it = 1,maxit
    efermi = 0.5d0*(e0 + e1) ! bisection
    chg = 0.d0
    DO ik = 1,nkpt
      DO ist = 1,nstsv
        e = evalsv(ist,ik)
        IF(e < e0min) THEN 
          occsv(ist,ik) = 0.d0
        else
          x = (efermi - e)*t1
          occsv(ist,ik) = occmax*stheta(stype,x)
          chg = chg + wkpt(ik)*occsv(ist,ik)
        ENDIF 
      ENDDO 
    ENDDO 
    IF(chg < chgval) THEN 
      e0 = efermi
    else
      e1 = efermi
    ENDIF 
    IF((e1-e0) < 1.d-12) goto 10 ! successful
  ENDDO 
  !
  WRITE(*,*)
  WRITE(*,'("Warning(occupy): could not find Fermi energy")')
  !
  10 continue
  ! find the density of states at the Fermi surface in units of
  ! states/Hartree/unit cell
  fermidos = 0.d0
  DO ik = 1,nkpt
    DO ist = 1,nstsv
      x = (evalsv(ist,ik) - efermi)*t1
      fermidos = fermidos + wkpt(ik)*sdelta(stype,x)*t1
    ENDDO 
    IF(abs(occsv(nstsv,ik)) > epsocc) THEN 
      WRITE(*,*)
      WRITE(*,'("Warning(occupy): not enough empty states for k-point ",I6)') ik
      WRITE(*,'(" and s.c. loop ",I5)') iscl
    ENDIF 
  ENDDO 
  fermidos = fermidos*occmax

  ! estimate the indirect band gap (FC)
  e0 = -1.d8
  e1 = 1.d8
  ikgap(1) = 1
  ikgap(2) = 1
  ! these loops are incorrectly ordered to fix a bug in versions 17 and 18 of the
  ! Intel compiler
  DO ist = 1,nstsv
    DO ik = 1,nkpt
      e = evalsv(ist,ik)
      IF(e < efermi) THEN 
        IF(e > e0) THEN 
          e0 = e
          ikgap(1) = ik
        ENDIF 
      else
        IF(e < e1) THEN 
          e1 = e
          ikgap(2) = ik
        ENDIF 
      ENDIF 
    ENDDO 
  ENDDO 
  bandgap(1) = e1 - e0

  ! write band gap to test file
  WRITE(*,*) 'Estimated indirect band gap: ', bandgap(1)

  ! estimate the direct band gap
  e = 1.d8
  ikgap(3) = 1
  DO ik = 1,nkpt
    e0 = -1.d8
    e1 = 1.d8
    DO ist = 1,nstsv
      t1 = evalsv(ist,ik)
      IF(t1 <= efermi) THEN 
        IF(t1 > e0) e0 = t1
      else
        IF(t1 < e1) e1 = t1
      ENDIF 
    ENDDO 
    t1 = e1 - e0
    IF(t1 < e) THEN 
      e = t1
      ikgap(3) = ik
    ENDIF 
  ENDDO 
  bandgap(2) = e

  WRITE(*,*) 'Estimated direct band gap: ', bandgap(1)

  RETURN 
END SUBROUTINE 
