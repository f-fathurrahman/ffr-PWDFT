
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvsig
! !INTERFACE:
subroutine genvsig
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the Kohn-Sham effective potential in the
!   interstitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
! multiply potential by characteristic function in real-space
zfft(:)=vsir(:)*cfunir(:)
! Fourier transform to G-space
call zfftifc(3,ngridg,-1,zfft)
! store in global array
vsig(1:ngvec)=zfft(igfft(1:ngvec))
deallocate(zfft)
return
end subroutine
!EOC

