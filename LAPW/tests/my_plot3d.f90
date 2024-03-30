! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} in the parallelepiped defined by the corner vertices in the
!   global array {\tt vclp3d}. See routine {\tt rfarray}.
SUBROUTINE my_plot3d(fnum, nf, rfmt, rfir)
  USE m_atoms, ONLY: natmtot
  USE m_gvectors, ONLY: ngtot
  USE m_lattice, ONLY: avec
  USE m_muffin_tins, ONLY: npmtmax
  USE m_plotting, ONLY: np3d, vclp3d
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: fnum,nf
  REAL(8), intent(in) :: rfmt(npmtmax,natmtot,nf), rfir(ngtot,nf)
  ! local variables
  INTEGER :: np,jf,ip
  REAL(8) :: v1(3)
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: vpl(:,:), fp(:,:)
  
  IF((nf < 1) .or. (nf > 4)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(plot3d): invalid number of functions : ",I8)') nf
    WRITE(*,*)
    STOP
  ENDIF 

  ! total number of plot points
  np = np3d(1)*np3d(2)*np3d(3)
  write(*,*) 'Total number of plot points = ', np

  ! allocate local arrays
  ALLOCATE(vpl(3,np),fp(np,nf))

  ! generate the 3D plotting points
  CALL my_plotpt3d(vclp3d, np3d, vpl)
  
  ! evaluate the functions at the grid points
  DO jf = 1,nf ! loop over number of functions
    CALL my_rfplot(np, vpl, rfmt(:,:,jf), rfir(:,jf), fp(:,jf))
  ENDDO 
  
  ! write functions to file
  WRITE(fnum,'(3I6," : grid size")') np3d(:)
  DO ip = 1,np
    CALL r3mv(avec, vpl(:,ip), v1)
    WRITE(fnum,'(7G18.10)') v1(:),(fp(ip,jf),jf=1,nf)
  ENDDO 
  DEALLOCATE(vpl,fp)
  RETURN 
END SUBROUTINE 