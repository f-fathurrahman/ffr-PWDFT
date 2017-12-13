! This module is a simplified version of ps_hgh.F90 in Octopus code

MODULE m_Ps_HGH

  USE m_constants, ONLY : PI

  IMPLICIT NONE

  TYPE Ps_HGH_Params_T
    CHARACTER(5) :: atom_name
    INTEGER :: zval
    REAL(8) :: rlocal
    REAL(8) :: rc(0:3)
    REAL(8) :: c(1:4)
    REAL(8) :: h(0:3, 1:3, 1:3)
    REAL(8) :: k(0:3, 1:3, 1:3)
    INTEGER :: lmax
    INTEGER :: Nproj_l(0:3)  ! number of projectors for each AM
    REAL(8) :: rcut_NL(0:3)
  END TYPE


CONTAINS


  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_G_long(p, g) RESULT(Vloc)
  !----------------------------------------------------------------------------
    type(Ps_HGH_Params_T), INTENT(in) :: p
    REAL(8), INTENT(in) :: g
    REAL(8) :: Vloc
    REAL(8), PARAMETER :: SMALL = 1.d-8

    REAL(8) :: g1, g2

    IF( g > SMALL ) THEN

      g1 = g*p%rlocal
      g2 = g1*g1

      Vloc = -(4.d0*PI*p%zval/g**2) * exp( -g2/2.d0)
    ELSE

      Vloc = 2.d0*PI * p%rlocal**2 * p%zval
    ENDIF

  END FUNCTION


  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_G(p, g) RESULT(Vloc)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T), INTENT(in) :: p
    REAL(8), INTENT(in) :: g
    REAL(8) :: Vloc
    REAL(8), PARAMETER :: SMALL = 1.d-8

    REAL(8) :: g1, g2, g4, g6

    IF( g > SMALL ) THEN

      g1 = g*p%rlocal
      g2 = g1*g1
      g4 = g2*g2
      g6 = g4*g2

      Vloc = -(4.d0*PI*p%zval/g**2) * exp( -g2/2.d0) + &
              sqrt(8.d0*PI**3) * p%rlocal**3 * exp( -g2/2.d0) * &
              ( p%c(1) + p%c(2)*(3.d0 - g2) + p%c(3)*(15.d0 - 10.d0*g2 + g4) + &
                p%c(4)*(105.d0 - 105.d0*g2 + 21.d0*g4 - g6) )
    ELSE

      Vloc = 2.d0*PI * p%rlocal**2 * p%zval + (2.d0*PI)**(1.5d0) * p%rlocal**3 * &
             ( p%c(1) + 3.d0*p%c(2) + 15.d0*p%c(3) + 105.d0*p%c(4) )

    ENDIF

  END FUNCTION


  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_proj_R_v2( psp, l, i, r ) RESULT( pil )
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    INTEGER :: l, i  ! l can be zero, i varies from 1 to 3
    REAL(8) :: r
    REAL(8) :: pil
    !
    REAL(8) :: num, denum
    !
    num = sqrt(2.d0)*r**(l + 2*(i-1)) * exp( -(r/psp%rc(l))**2 )
    denum = (psp%rc(l))**( l + (4*i-1)/2 ) * gamma( dble(l) + (4*i-1)/2.d0 )
    pil = num/denum
  END FUNCTION


  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_proj_R(p, l, i, r) RESULT(f_prj)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T), INTENT(in) :: p
    REAL(8), INTENT(in) :: r
    INTEGER, INTENT(in) :: i, l
    REAL(8) :: f_prj

    REAL(8) :: x, y, rr

    x = l + real(4*i-1, kind=8)/2.d0
    y = gamma(x)  ! use intrinsic function ?
    x = sqrt(y)
    IF(l==0 .and. i==1) THEN
      rr = 1.d0
    ELSE
      rr = r ** (l + 2*(i-1))
    ENDIF

    f_prj = sqrt(2.d0) * rr * exp(-r**2/(2.d0*p%rc(l)**2)) / &
            ( p%rc(l)**(l + real(4*i-1, kind=8)/2.d0) * x )

  END FUNCTION


  ! evaluate 4*pi*r2*Vloc(r)
  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_4pi_r2( psp, r ) RESULT(Vloc)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    REAL(8) :: r
    !
    INTEGER :: i
    REAL(8) :: rrloc ! r/rlocal
    REAL(8) :: term1, Vloc

    rrloc = r/psp%rlocal
    term1 = psp%c(1)
    DO i = 2, 4
      term1 = term1 + psp%c(i)*rrloc**(2*(i-1))
    ENDDO
    Vloc = -4.d0*PI*r*psp%zval * erf( rrloc/sqrt(2.d0) ) + exp(-0.5d0*rrloc**2)*term1 * 4.d0*pi*r**2
  END FUNCTION

  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_R( psp, r ) RESULT(Vloc)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    REAL(8) :: r
    !
    INTEGER :: i
    REAL(8) :: rrloc ! r/rlocal
    REAL(8) :: term1, Vloc
    REAL(8), PARAMETER :: SMALL = 1d-10

    IF ( r <= SMALL ) THEN
      rrloc = SMALL/psp%rlocal
      term1 = psp%c(1)
      DO i = 2, 4
        term1 = term1 + psp%c(i)*rrloc**(2*(i-1))
      ENDDO
      Vloc = -psp%zval/SMALL * erf( rrloc/sqrt(2.d0) ) + exp(-0.5d0*rrloc**2)*term1
    ELSE
      rrloc = r/psp%rlocal
      term1 = psp%c(1)
      DO i = 2, 4
        term1 = term1 + psp%c(i)*rrloc**(2*(i-1))
      ENDDO
      Vloc = -psp%zval/r * erf( rrloc/sqrt(2.d0) ) + exp(-0.5d0*rrloc**2)*term1
    ENDIF
  END FUNCTION

  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_R_long( psp, r ) RESULT(Vloc)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    REAL(8) :: r
    !
    REAL(8) :: rrloc ! r/rlocal
    REAL(8) :: Vloc
    REAL(8), PARAMETER :: SMALL = 1d-10

    IF ( r <= SMALL ) THEN
      rrloc = SMALL/psp%rlocal
      Vloc = -psp%zval/SMALL * erf( rrloc/sqrt(2.d0) )
    ELSE
      rrloc = r/psp%rlocal
      Vloc = -psp%zval/r * erf( rrloc/sqrt(2.d0) )
    ENDIF
  END FUNCTION


  !----------------------------------------------------------------------------
  FUNCTION hgh_eval_Vloc_R_short( psp, r ) RESULT(Vloc)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    REAL(8) :: r
    !
    INTEGER :: i
    REAL(8) :: rrloc ! r/rlocal
    REAL(8) :: term1, Vloc

    rrloc = r/psp%rlocal
    term1 = psp%c(1)
    DO i = 2, 4
      term1 = term1 + psp%c(i)*rrloc**(2*(i-1))
    ENDDO
    Vloc = exp(-0.5d0*rrloc**2)*term1
  END FUNCTION


  ! format used by CP2K program (slightly modified)
  !-----------------------------------------------
  SUBROUTINE init_Ps_HGH_Params(psp, filename)
  !-----------------------------------------------
    TYPE(Ps_HGH_Params_T), INTENT(INOUT) :: psp
    CHARACTER(len=*), INTENT(IN) :: filename
    !
    INTEGER :: iunit, ii, i, j, k
    LOGICAL :: found
    INTEGER :: n_s, n_p, n_d, n_f, n_c_local

    ! Set initially everything to zero.
    psp%c(1:4)  = 0.d0
    psp%rlocal  = 0.d0
    psp%rc(0:3) = 0.d0
    psp%h(0:3,1:3,1:3) = 0.d0
    psp%k(0:3,1:3,1:3) = 0.d0
    psp%Nproj_l(0:3) = 0

    INQUIRE(file=filename, exist=found)
    IF(.NOT.found) THEN
      WRITE(*,*) 'ERROR: Pseudopotential file: ', trim(filename), ' is not found.'
      STOP
    ENDIF

    iunit = 303  ! arbirary
    OPEN(unit=iunit,file=filename, action='read', form='formatted', status='old')

    READ(iunit,*) psp%atom_name
    READ(iunit,*) n_s, n_p, n_d, n_f
    psp%zval = n_s + n_p + n_d + n_f
    
    READ(iunit,*) psp%rlocal, n_c_local
    !WRITE(*,*) 'n_c_local = ', n_c_local
    READ(iunit,*) psp%c(1:n_c_local)

    READ(iunit,*) psp%lmax
    psp%lmax = psp%lmax - 1
    !WRITE(*,*) 'psp%lmax = ', psp%lmax

    DO i = 0,psp%lmax
      READ(iunit,*) psp%rc(i), psp%Nproj_l(i)
      !WRITE(*,*) psp%rc(i), psp%Nproj_l(i)
      DO ii = 1,psp%Nproj_l(i)
        READ(iunit,*) psp%h(i,ii,ii:psp%Nproj_l(i))
      ENDDO 
    ENDDO 
    
    CLOSE(iunit)

    ! Parameters are symmetric.
    DO k = 0, 3
      DO i = 1, 3
        DO j = i + 1, 3
          psp%h(k, j, i) = psp%h(k, i, j)
        ENDDO
      ENDDO
    ENDDO

    CALL find_NL_cutoff(psp)

  END SUBROUTINE 



  !----------------------------------------------------------------------------
  SUBROUTINE init_Ps_HGH_Params_old(psp, filename)
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T), INTENT(INOUT) :: psp
    CHARACTER(len=*), INTENT(IN)  :: filename

    INTEGER :: iunit, i, l
    LOGICAL :: found

    INQUIRE(file=filename, exist=found)
    IF(.NOT.found) THEN
      WRITE(*,*) 'ERROR: Pseudopotential file: ', trim(filename), ' is not found.'
      STOP
    ENDIF

    iunit = 303  ! just saying
    OPEN(unit=iunit,file=filename, action='read', form='formatted', status='old')
    i = load_params(iunit, psp)
    IF(i /= 0) THEN
      WRITE(*,*) 'ERROR: Reading psfile: ', trim(filename), '.'
      STOP
    ENDIF
    CLOSE(iunit)

    ! Find Nproj_l based on values of h
    DO l = 0,3
      DO i = 1,3
        IF( abs( psp%h(l,i,i) ) > 0.d0 ) psp%Nproj_l(l) = psp%Nproj_l(l) + 1
      ENDDO
      !WRITE(*,*) l, psp%Nproj_l(l)
    ENDDO

    ! Find Nproj_l based on values of rc if Nproj_l == 0
    !DO l = 0,psp%lmax
    !  IF( psp%Nproj_l(l) == 0 ) THEN
    !    IF( abs(psp%rc(l)) > 0.d0 ) psp%Nproj_l(l) = 1
    !  ENDIF
    !ENDDO

    ! Finds out psp%lmax. The most special cases are H, He, Li_sc and Be_sc, where psp%lmax = -1.
    !psp%lmax = 0
    !DO WHILE(psp%rc(psp%lmax) > 0.01d0)
    !  psp%lmax = psp%lmax + 1
    !  IF(psp%lmax > 3) EXIT
    !ENDDO
    !psp%lmax = psp%lmax - 1

    psp%lmax = -1
    DO l = 0,3
      IF( psp%Nproj_l(l) > 0 ) psp%lmax = psp%lmax + 1
    ENDDO

    CALL find_NL_cutoff(psp)

  END SUBROUTINE


  !----------------------------------------------------------------------------
  FUNCTION load_params(iunit, params)
  !----------------------------------------------------------------------------
    USE m_xc, ONLY : XC_NAME
    INTEGER, INTENT(IN)  :: iunit ! where to read from
    TYPE(Ps_HGH_Params_T), INTENT(out) :: params ! should INOUT instead?
    INTEGER :: load_params ! 0 if success, 1 otherwise.
    CHARACTER(80) :: line

    INTEGER :: i, ios, j, k

    ! Set initially everything to zero.
    params%c(1:4) = 0.d0
    params%rlocal = 0.d0
    params%rc = 0.d0
    params%h = 0.d0
    params%k = 0.d0
    params%Nproj_l(:) = 0

    ! get valence configuration
    READ(iunit,'(a)') line
    !WRITE(*,*) 'Valconf = ', line

    ! Reads the file in a hopefully smart way
    ios = 1
    j = 5
    READ(iunit,'(a)') line
    DO WHILE((ios /= 0) .AND. (j > 0))
      j = j - 1
      READ(line, *, iostat=ios) params%atom_name, params%zval, params%rlocal, params%c(1:j)
    ENDDO
    IF (j<1) READ(line, *, iostat=ios) params%atom_name, params%zval, params%rlocal
    IF ( ios /= 0 ) then
      load_params = 1
      RETURN
    ENDIF

    READ(iunit,'(a)', iostat = ios) line
    IF(ios /= 0) THEN
      load_params = 0
      RETURN
    ENDIF
    ios = 1
    j = 4
    DO WHILE((ios /= 0) .AND. (j > 0))
      j = j - 1
      READ(line, *, iostat=ios) params%rc(0), (params%h(0, i, i), i = 1, j)
    ENDDO
    IF(j < 0) THEN
      load_params = 2
      RETURN
    ENDIF

    kloop: DO k = 1, 3
      READ(iunit, '(a)', iostat=ios) line
      IF(ios /= 0) EXIT kloop
      ios = 1
      j = 4
      DO WHILE((ios /= 0) .AND. (j > 0))
        j = j - 1
        READ(line, *, iostat=ios) params%rc(k), (params%h(k, i, i), i = 1, j)
      ENDDO
      IF(params%rc(k) == 0.d0) EXIT kloop
      READ(iunit, '(a)') line
      ios = 1
      j = 4
      DO WHILE((ios /= 0) .AND. (j>0))
        j = j - 1
        READ(line, *, iostat=ios) (params%k(k, i, i), i = 1, 3)
      ENDDO
    ENDDO kloop

    ! Fill in the rest of the parameter matrices...
    params%h(0, 1, 2) = -0.5d0      * sqrt(3.d0/5.d0)       * params%h(0, 2, 2)
    params%h(0, 1, 3) =  0.5d0      * sqrt(5.d0/21.d0)      * params%h(0, 3, 3)
    params%h(0, 2, 3) = -0.5d0      * sqrt(100.d0/63.d0)    * params%h(0, 3, 3)
    params%h(1, 1, 2) = -0.5d0      * sqrt(5.d0/7.d0)       * params%h(1, 2, 2)
    params%h(1, 1, 3) =  1.d0/6.d0  * sqrt(35.d0/11.d0)     * params%h(1, 3, 3)
    params%h(1, 2, 3) = -1.d0/6.d0  * (14.d0 / sqrt(11.d0)) * params%h(1, 3, 3)
    params%h(2, 1, 2) = -0.5d0      * sqrt(7.d0/9.d0)       * params%h(2, 2, 2)
    params%h(2, 1, 3) =  0.5d0      * sqrt(63.d0/143.d0)    * params%h(2, 3, 3)
    params%h(2, 2, 3) = -0.5d0      * (18.d0/sqrt(143.0d0)) * params%h(2, 3, 3)

    params%k(0, 1, 2) = -0.5d0      * sqrt(3.d0/5.d0)       * params%k(0, 2, 2)
    params%k(0, 1, 3) =  0.5d0      * sqrt(5.d0/21.d0)      * params%k(0, 3, 3)
    params%k(0, 2, 3) = -0.5d0      * sqrt(100.d0/63.d0)    * params%k(0, 3, 3)
    params%k(1, 1, 2) = -0.5d0      * sqrt(5.d0/7.d0)       * params%k(1, 2, 2)
    params%k(1, 1, 3) =  1.d0/6.d0  * sqrt(35.d0/11.d0)     * params%k(1, 3, 3)
    params%k(1, 2, 3) = -1.d0/6.d0  * (14.d0/sqrt(11.d0))   * params%k(1, 3, 3)
    params%k(2, 1, 2) = -0.5d0      * sqrt(7.d0/9.d0)       * params%k(2, 2, 2)
    params%k(2, 1, 3) =  0.5d0      * sqrt(63.0/143.0)      * params%k(2, 3, 3)
    params%k(2, 2, 3) = -0.5d0      * (18.d0/sqrt(143.d0))  * params%k(2, 3, 3)

    ! special case for PBE ????
    IF( XC_NAME == 'PBE' ) THEN 
      params%h(0, 1, 2) = 2.d0*params%h(0, 1, 2)
      params%h(0, 1, 3) = 2.d0*params%h(0, 1, 3)
      params%h(0, 2, 3) = 2.d0*params%h(0, 2, 3)
      params%h(1, 1, 2) = 2.d0*params%h(1, 1, 2)
      params%h(1, 1, 3) = 2.d0*params%h(1, 1, 3)
      params%h(1, 2, 3) = 2.d0*params%h(1, 2, 3)
      params%h(2, 1, 2) = 2.d0*params%h(2, 1, 2)
      params%h(2, 1, 3) = 2.d0*params%h(2, 1, 3)
      params%h(2, 2, 3) = 2.d0*params%h(2, 2, 3)
    ENDIF 

    ! Parameters are symmetric.
    DO k = 0, 3
      DO i = 1, 3
        DO j = i + 1, 3
          params%h(k, j, i) = params%h(k, i, j)
          params%k(k, j, i) = params%k(k, i, j)
        ENDDO
      ENDDO
    ENDDO

    load_params = 0
  END FUNCTION load_params


  ! To mimic functionality of get_cutoff_radii
  !----------------------------------------------------------------------------
  SUBROUTINE find_NL_cutoff( psp )
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: psp
    !
    INTEGER :: i, iprj, l
    REAL(8) :: r, rcut_tmp
    !
    INTEGER, PARAMETER :: NRMAX = 2000
    REAL(8), PARAMETER :: dr = 0.01d0
    REAL(8), PARAMETER :: SMALL = 1.d-5
    !
    REAL(8) :: f_prj

    psp%rcut_NL(:) = 0.d0
    rcut_tmp = 0.d0

    DO l = 0, psp%lmax
      DO iprj = 1, psp%Nproj_l(l)
        DO i = 1,NRMAX
          r = psp%rc(l) + i*dr
          f_prj = hgh_eval_proj_R( psp, l, iprj, r )
          IF( f_prj < SMALL ) THEN
            rcut_tmp = r
            !WRITE(*,*) 'l, iprj, i, r', l, iprj, i, rcut_tmp
            EXIT
          ENDIF
        ENDDO
        psp%rcut_NL(l) = max( psp%rcut_NL(l), rcut_tmp )
      ENDDO
    ENDDO

  END SUBROUTINE


  !----------------------------------------------------------------------------
  SUBROUTINE info_Ps_HGH_Params( ps )
  !----------------------------------------------------------------------------
    TYPE(Ps_HGH_Params_T) :: ps
    INTEGER :: i, j, l

    WRITE(*,'(4x,2A)') 'atom_name = ', trim(ps%atom_name)
    WRITE(*,'(4x,A,I5)') 'zval = ', ps%zval
    WRITE(*,'(4x,A,I5)') 'lmax = ', ps%lmax
    WRITE(*,*)
    WRITE(*,'(4x,A)') 'Local pseudopotential parameters:'
    WRITE(*,*)
    WRITE(*,'(4x,A,F18.10)') 'rlocal = ', ps%rlocal
    !
    WRITE(*,*)
    WRITE(*,'(4x,A)') 'c (for local potential ) = '
    DO i = 1, 4
      WRITE(*,'(4x,I5,F18.10)') i, ps%c(i)
    ENDDO

    IF( ps%lmax >= 0 ) THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') 'Nonlocal pseudopotential parameters:'
      DO l = 0, ps%lmax
        WRITE(*,*)
        WRITE(*,'(4x,A,I4)') 'Matrix h for l = ', l
        WRITE(*,'(4x,A,I4)') 'Nproj_l        = ', ps%Nproj_l(l)
        WRITE(*,'(4x,A,F18.10)') 'rc             = ', ps%rc(l)
        WRITE(*,'(4x,A,F18.10)') 'rcut_NL        = ', ps%rcut_NL(l)
        DO i = 1, 3
          WRITE(*,*)
          DO j = 1, 3
            WRITE(*,'(F18.10)',advance='no') ps%h(l,i,j)
          ENDDO
        ENDDO
        WRITE(*,*)
      ENDDO

    ENDIF
  END SUBROUTINE


END MODULE
