module m_sht

!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! trotsht is .true. if the spherical cover used for the SHT is to be rotated
logical trotsht
data trotsht / .false. /
! spherical cover rotation matrix
real(8) rotsht(3,3)

! real backward SHT matrix for lmaxi
real(8), allocatable :: rbshti(:,:)

! real forward SHT matrix for lmaxi
real(8), allocatable :: rfshti(:,:)

! real backward SHT matrix for lmaxo
real(8), allocatable :: rbshto(:,:)

! real forward SHT matrix for lmaxo
real(8), allocatable :: rfshto(:,:)

! complex backward SHT matrix for lmaxi
complex(8), allocatable :: zbshti(:,:)

! complex forward SHT matrix for lmaxi
complex(8), allocatable :: zfshti(:,:)

! complex backward SHT matrix for lmaxo
complex(8), allocatable :: zbshto(:,:)

! complex forward SHT matrix for lmaxo
complex(8), allocatable :: zfshto(:,:)

end module
