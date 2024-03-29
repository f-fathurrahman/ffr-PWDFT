! !DESCRIPTION:
!   Generates the smooth characteristic function. This is the function which is
!   0 within the muffin-tins and 1 in the intersitial region and is constructed
!   from radial step function form factors with $G<G_{\rm max}$. The form
!   factors are given by
!   $$ \tilde{\Theta}_i(G)=\begin{cases}
!    \frac{4\pi R_i^3}{3 \Omega} & G=0 \\
!    \frac{4\pi R_i^3}{\Omega}\frac{j_1(GR_i)}{GR_i} & 0<G\le G_{\rm max} \\
!    0 & G>G_{\rm max}\end{cases} $$
!   where $R_i$ is the muffin-tin radius of the $i$th species and $\Omega$ is
!   the unit cell volume. Therefore the characteristic function in $G$-space is
!   $$ \tilde{\Theta}({\bf G})=\delta_{G,0}-\sum_{ij}\exp(-i{\bf G}\cdot
!    {\bf r}_{ij})\tilde{\Theta}_i(G), $$
!   where ${\bf r}_{ij}$ is the position of the $j$th atom of the $i$th species.
SUBROUTINE gencfun()
  USE m_atoms, ONLY: natoms, nspecies, atposc
  USE m_gvectors, ONLY: cfunig, cfunir, ffacg, ngtot, ngridg, igfft, vgc
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ig
  REAL(8) :: v1,v2,v3,t1
  COMPLEX(8) :: z1
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: zfft(:)
  ! allocate global characteristic function arrays
  IF(allocated(cfunig)) DEALLOCATE(cfunig)
  ALLOCATE(cfunig(ngtot))
  IF(allocated(cfunir)) DEALLOCATE(cfunir)
  ALLOCATE(cfunir(ngtot))
  cfunig(1)=1.d0
  cfunig(2:)=0.d0
  ! begin loop over species
  DO is=1,nspecies
    ! loop over atoms
    DO ia=1,natoms(is)
      v1=atposc(1,ia,is); v2=atposc(2,ia,is); v3=atposc(3,ia,is)
      DO ig=1,ngtot
        ! structure factor
        t1=vgc(1,ig)*v1+vgc(2,ig)*v2+vgc(3,ig)*v3
        z1=cmplx(cos(t1),-sin(t1),8)
        ! add to characteristic function in G-space
        cfunig(ig)=cfunig(ig)-ffacg(ig,is)*z1
      ENDDO 
    ENDDO 
  ENDDO 
  ALLOCATE(zfft(ngtot))
  zfft(igfft(:))=cfunig(:)
  ! Fourier transform to real-space
  CALL zfftifc(3,ngridg,1,zfft)
  cfunir(:)=dble(zfft(:))
  DEALLOCATE(zfft)
  RETURN 
END SUBROUTINE 
