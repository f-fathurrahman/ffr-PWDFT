! Calculate structure factor
SUBROUTINE init_strfact()
 
  USE m_atoms, ONLY : Na => Natoms, &
                      Xpos => AtomicCoords, &
                      Nspecies, &
                      atm2species, &
                      strf => StructureFactor
  USE m_PWGrid, ONLY : Ng, Gv
  IMPLICIT NONE
  !
  INTEGER :: ia, isp, ig
  REAL(8) :: GX

  WRITE(*,*)
  WRITE(*,*) 'Calculating structure factor'

  ALLOCATE( strf(Ng,Nspecies) )

  strf(:,:) = cmplx(0.d0,0.d0,kind=8)
  DO ia = 1,Na
    isp = atm2species(ia)
    DO ig = 1,Ng
      GX = Xpos(1,ia)*Gv(1,ig) + Xpos(2,ia)*Gv(2,ig) + Xpos(3,ia)*Gv(3,ig)
      strf(ig,isp) = strf(ig,isp) + cmplx( cos(GX), -sin(GX), kind=8 )
    ENDDO 
  ENDDO 

END SUBROUTINE 

