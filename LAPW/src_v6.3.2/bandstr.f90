
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:
subroutine bandstr
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Produces a band structure along the path in reciprocal space which connects
!   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
!   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
!   with the Fermi energy set to zero. If required, band structures are plotted
!   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
!   which include the band characters for each $l$ component of that atom in
!   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
!   Vertex location lines are written to {\tt BANDLINES.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn,is,ia,ias
integer lmax,lmmax,l,m,lm,iv,nthd
real(8) emin,emax,sum
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:),e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:),apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! initialise universal variables
call init0
call init1
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
! maximum angular momentum for band character
lmax=min(3,lmaxo)
lmmax=(lmax+1)**2
if (task.eq.21) then
  allocate(bc(0:lmax,natmtot,nstsv,nkpt))
else if (task.eq.22) then
  allocate(bc(lmmax,natmtot,nstsv,nkpt))
else
  allocate(bc(nspinor,natmtot,nstsv,nkpt))
end if
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr
emin=1.d5
emax=-1.d5
! begin parallel loop over k-points
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP PRIVATE(dmat,apwalm,ist,ispn) &
!$OMP PRIVATE(ias,l,m,lm,sum) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
if (task.ge.21) then
  allocate(dmat(lmmax,nspinor,lmmax,nspinor,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
end if
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL(bandstr_1)
  write(*,'("Info(bandstr): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(bandstr_1)
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ik)=evalsv(ist,ik)-efermi
!$OMP CRITICAL(bandstr_2)
    emin=min(emin,e(ist,ik))
    emax=max(emax,e(ist,ik))
!$OMP END CRITICAL(bandstr_2)
  end do
! compute the band characters if required
  if (task.ge.21) then
! find the matching coefficients
    do ispn=1,nspnfv
      call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    end do
! average band character over spin and m for all atoms
    do ias=1,natmtot
! generate the diagonal of the density matrix
      call gendmatk(.true.,.true.,0,lmax,ias,ngk(:,ik),apwalm,evecfv,evecsv, &
       lmmax,dmat)
      do ist=1,nstsv
        if (task.eq.21) then
! l character of band
          lm=0
          do l=0,lmax
            sum=0.d0
            do m=-l,l
              lm=lm+1
              do ispn=1,nspinor
                sum=sum+dble(dmat(lm,ispn,lm,ispn,ist))
              end do
            end do
            bc(l,ias,ist,ik)=real(sum)
          end do
        else if (task.eq.22) then
! (l,m) character of band
          lm=0
          do l=0,lmax
            do m=-l,l
              lm=lm+1
              sum=0.d0
              do ispn=1,nspinor
                sum=sum+dble(dmat(lm,ispn,lm,ispn,ist))
              end do
              bc(lm,ias,ist,ik)=real(sum)
            end do
          end do
        else
! spin character of band
          do ispn=1,nspinor
            sum=0.d0
            lm=0
            do l=0,lmax
              do m=-l,l
                lm=lm+1
                sum=sum+dble(dmat(lm,ispn,lm,ispn,ist))
              end do
            end do
            bc(ispn,ias,ist,ik)=real(sum)
          end do
        end if
      end do
    end do
  end if
! end loop over k-points
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)
if (task.ge.21) deallocate(dmat,apwalm)
!$OMP END PARALLEL
call freethd(nthd)
emax=emax+(emax-emin)*0.5d0
emin=emin-(emax-emin)*0.5d0
! output the band structure
if (task.eq.20) then
  open(50,file='BAND.OUT',form='FORMATTED')
  do ist=1,nstsv
    do ik=1,nkpt
      write(50,'(2G18.10)') dpp1d(ik),e(ist,ik)
    end do
    write(50,'("     ")')
  end do
  close(50)
  write(*,*)
  write(*,'("Info(bandstr):")')
  write(*,'(" band structure plot written to BAND.OUT")')
else
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(fname,'("BAND_S",I2.2,"_A",I4.4,".OUT")') is,ia
      open(50,file=trim(fname),form='FORMATTED')
      do ist=1,nstsv
        do ik=1,nkpt
          if (task.eq.21) then
! sum band character over l to find total atomic character
            sum=0.d0
            do l=0,lmax
              sum=sum+bc(l,ias,ist,ik)
            end do
            write(50,'(2G18.10,5F12.6)') dpp1d(ik),e(ist,ik),sum, &
             (bc(l,ias,ist,ik),l=0,lmax)
          else if (task.eq.22) then
            write(50,'(2G18.10,16F12.6)') dpp1d(ik),e(ist,ik), &
             (bc(lm,ias,ist,ik),lm=1,lmmax)
          else
            write(50,'(2G18.10,2F12.6)') dpp1d(ik),e(ist,ik), &
             (bc(ispn,ias,ist,ik),ispn=1,nspinor)
          end if
        end do
        write(50,'("     ")')
      end do
      close(50)
    end do
  end do
  write(*,*)
  write(*,'("Info(bandstr):")')
  write(*,'(" Band structure plot written to BAND_Sss_Aaaaa.OUT")')
  write(*,'(" for all species and atoms")')
  write(*,*)
  write(*,'(" Columns in the file are :")')
  if (task.eq.21) then
    write(*,'("  distance, eigenvalue, total atomic character, l character &
     &(l = 0...",I1,")")') lmax
  else if (task.eq.22) then
    write(*,'("  distance, eigenvalue, (l,m) character &
     &(l = 0...",I1,", m = -l...l)")') lmax
  else
    write(*,'("  distance, eigenvalue, spin-up and spin-down characters")')
  end if
end if
write(*,*)
write(*,'(" Fermi energy is at zero in plot")')
! output the vertex location lines
open(50,file='BANDLINES.OUT',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'(" Vertex location lines written to BANDLINES.OUT")')
deallocate(e)
if (task.ge.21) deallocate(bc)
return
end subroutine
!EOC

