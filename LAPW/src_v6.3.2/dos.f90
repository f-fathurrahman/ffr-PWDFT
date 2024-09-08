
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dos
! !INTERFACE:
subroutine dos(fext,tocc,occsvp)
! !USES:
use modmain
use modomp
use modtest
! !INPUT/OUTPUT PARAMETERS:
!   fext   : filename extension (in,character(*))
!   tocc   : .true. if just the occupied orbitals should contribute to the DOS
!            (in,logical)
!   occsvp : occupation numbers of second-variational orbitals
!            (in,real(nstsv,nkpt))
! !DESCRIPTION:
!   Produces a total and partial density of states (DOS) for plotting. The total
!   DOS is written to the file {\tt TDOS.OUT} while the partial DOS is written
!   to the file {\tt PDOS\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species
!   {\tt ss}. In the case of the partial DOS, each symmetrised
!   $(l,m)$-projection is written consecutively and separated by blank lines.
!   If the global variable {\tt lmirep} is {\tt .true.}, then the density matrix
!   from which the $(l,m)$-projections are obtained is first rotated into a
!   irreducible representation basis, i.e. one that block diagonalises all the
!   site symmetry matrices in the $Y_{lm}$ basis. Eigenvalues of a quasi-random
!   matrix in the $Y_{lm}$ basis which has been symmetrised with the site
!   symmetries are written to {\tt ELMIREP.OUT}. This allows for identification
!   of the irreducible representations of the site symmetries, for example $e_g$
!   or $t_{2g}$, by the degeneracies of the eigenvalues. In the plot, spin-up is
!   made positive and spin-down negative. See the routines {\tt gendmatk} and
!   {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Parallelised and included sum over m, November 2009 (F. Cricchio)
!EOP
!BOC
implicit none
! arguments
character(*), intent(in) :: fext
logical, intent(in) :: tocc
real(8), intent(in) :: occsvp(nstsv,nkpt)
! local variables
logical tsqaz
integer nsk(3),ik,jk,ist,iw,ld
integer nsd,ispn,jspn,is,ia,ias
integer lmmax,l0,l1,l,m,lm,nthd
real(8) dw,th,sps(2),vl(3),vc(3)
real(8) v1(3),v2(3),v3(3),t1
complex(8) su2(2,2),b(2,2),c(2,2)
character(256) fname
! low precision for band/spin character array saves memory
real(4), allocatable :: bc(:,:,:,:,:),sc(:,:,:)
real(8), allocatable :: w(:),e(:,:,:),f(:,:)
real(8), allocatable :: g(:),dt(:,:),dp(:,:,:)
real(8), allocatable :: elm(:,:)
complex(8), allocatable :: ulm(:,:,:),a(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:),sdmat(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
lmmax=(lmaxdos+1)**2
ld=lmmax*nspinor
if (dosssum) then
  nsd=1
else
  nsd=nspinor
end if
if (dosmsum) then
  l0=0; l1=lmaxdos
else
  l0=1; l1=lmmax
end if
allocate(bc(lmmax,nspinor,natmtot,nstsv,nkptnr))
allocate(sc(nspinor,nstsv,nkptnr))
! generate unitary matrices which convert the (l,m) basis into the irreducible
! representation basis of the symmetry group at each atomic site
if (lmirep) then
  allocate(elm(lmmax,natmtot))
  allocate(ulm(lmmax,lmmax,natmtot))
  call genlmirep(lmaxdos,lmmax,elm,ulm)
end if
! compute the SU(2) operator used for rotating the density matrix to the
! desired spin-quantisation axis
v1(:)=sqados(:)
t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
if (t1.le.epslat) then
  write(*,*)
  write(*,'("Error(dos): spin-quantisation axis (sqados) has zero length")')
  write(*,*)
  stop
end if
v1(:)=v1(:)/t1
if (v1(3).ge.1.d0-epslat) then
  tsqaz=.true.
else
  tsqaz=.false.
  v2(1:2)=0.d0
  v2(3)=1.d0
  call r3cross(v1,v2,v3)
! note that the spin-quantisation axis is rotated, so the density matrix should
! be rotated in the opposite direction
  th=-acos(v1(3))
  call axangsu2(v3,th,su2)
end if
! begin parallel loop over k-points
call holdthd(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,dmat,sdmat,a) &
!$OMP PRIVATE(jk,ispn,jspn,vl,vc) &
!$OMP PRIVATE(is,ia,ias,ist,lm,b,c,t1) &
!$OMP NUM_THREADS(nthd)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(dmat(lmmax,nspinor,lmmax,nspinor,nstsv))
allocate(sdmat(nspinor,nspinor,nstsv),a(lmmax,lmmax))
!$OMP DO
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! loop over first-variational spins
  do ispn=1,nspnfv
    vl(:)=vkl(:,ik)
    vc(:)=vkc(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (ispn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! find the matching coefficients
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file for non-reduced k-point
  call getevecfv('.OUT',0,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv('.OUT',0,vkl(:,ik),evecsv)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! generate the density matrix
      call gendmatk(.false.,.false.,0,lmaxdos,ias,ngk(:,ik),apwalm,evecfv, &
       evecsv,lmmax,dmat)
! convert (l,m) part to an irreducible representation if required
      if (lmirep) then
        do ist=1,nstsv
          do ispn=1,nspinor
            do jspn=1,nspinor
              call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulm(:,:,ias),lmmax, &
               dmat(:,ispn,1,jspn,ist),ld,zzero,a,lmmax)
              call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulm(:,:,ias), &
               lmmax,zzero,dmat(:,ispn,1,jspn,ist),ld)
            end do
          end do
        end do
      end if
! spin rotate the density matrices to desired spin-quantisation axis
      if (spinpol.and.(.not.tsqaz)) then
        do ist=1,nstsv
          do lm=1,lmmax
            b(:,:)=dmat(lm,:,lm,:,ist)
            call z2mm(su2,b,c)
            call z2mmct(c,su2,b)
            dmat(lm,:,lm,:,ist)=b(:,:)
          end do
        end do
      end if
! determine the band characters from the density matrix
      do ist=1,nstsv
        do ispn=1,nspinor
          do lm=1,lmmax
            t1=dble(dmat(lm,ispn,lm,ispn,ist))
            bc(lm,ispn,ias,ist,ik)=real(t1)
          end do
        end do
      end do
    end do
  end do
! compute the spin density matrices of the second-variational states
  call gensdmat(evecsv,sdmat)
! spin rotate the density matrices to desired spin-quantisation axis
  if (spinpol.and.(.not.tsqaz)) then
    do ist=1,nstsv
      call z2mm(su2,sdmat(:,:,ist),b)
      call z2mmct(b,su2,sdmat(:,:,ist))
    end do
  end if
  do ist=1,nstsv
    do ispn=1,nspinor
      t1=dble(sdmat(ispn,ispn,ist))
      sc(ispn,ist,ik)=real(t1)
    end do
  end do
end do
!$OMP END DO
deallocate(apwalm,evecfv,evecsv,dmat,sdmat,a)
!$OMP END PARALLEL
call freethd(nthd)
allocate(w(nwplot),e(nstsv,nkptnr,nspinor))
allocate(dt(nwplot,nsd),dp(nwplot,l0:l1,nsd))
! generate frequency grid
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  w(iw)=dw*dble(iw-1)+wplot(1)
end do
! number of subdivisions used for interpolation in the Brillouin zone
nsk(:)=max(ngrkf/ngridk(:),1)
! sign for spin in DOS
sps(1)=1.d0
sps(2)=-1.d0
!-------------------!
!     total DOS     !
!-------------------!
allocate(f(nstsv,nkptnr),g(nwplot))
dt(:,:)=0.d0
do ispn=1,nspinor
  do ik=1,nkptnr
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
    do ist=1,nstsv
! subtract the Fermi energy
      e(ist,ik,ispn)=evalsv(ist,jk)-efermi
! use diagonal of spin density matrix for weight
      f(ist,ik)=sc(ispn,ist,ik)
      if (tocc) then
        f(ist,ik)=f(ist,ik)*occsvp(ist,jk)
      else
        f(ist,ik)=f(ist,ik)*occmax
      end if
    end do
  end do
! integrate over the Brillouin zone
  call brzint(nswplot,ngridk,nsk,ivkiknr,nwplot,wplot,nstsv,nstsv,e(:,:,ispn), &
   f,g)
  if (dosssum) then
    dt(:,1)=dt(:,1)+g(:)
  else
    dt(:,ispn)=g(:)
  end if
end do
deallocate(f,g)
! output to file
open(50,file='TDOS'//trim(fext),form='FORMATTED')
do ispn=1,nsd
  do iw=1,nwplot
    write(50,'(2G18.10)') w(iw),dt(iw,ispn)*sps(ispn)
  end do
  write(50,'("     ")')
end do
close(50)
!---------------------!
!     partial DOS     !
!---------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    dp(:,:,:)=0.d0
    do ispn=1,nspinor
      call holdthd(lmmax,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(f,g,l,ik,jk,ist) &
!$OMP NUM_THREADS(nthd)
      allocate(f(nstsv,nkptnr),g(nwplot))
!$OMP DO
      do lm=1,lmmax
        l=idxil(lm)
        do ik=1,nkptnr
          jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
          do ist=1,nstsv
            f(ist,ik)=bc(lm,ispn,ias,ist,ik)
            if (tocc) then
              f(ist,ik)=f(ist,ik)*occsvp(ist,jk)
            else
              f(ist,ik)=f(ist,ik)*occmax
            end if
          end do
        end do
        call brzint(nswplot,ngridk,nsk,ivkiknr,nwplot,wplot,nstsv,nstsv, &
         e(:,:,ispn),f,g)
        if (dosmsum) then
          if (dosssum) then
            dp(:,l,1)=dp(:,l,1)+g(:)
          else
            dp(:,l,ispn)=dp(:,l,ispn)+g(:)
          end if
        else
          if (dosssum) then
            dp(:,lm,1)=dp(:,lm,1)+g(:)
          else
            dp(:,lm,ispn)=g(:)
          end if
        end if
! subtract from interstitial DOS
!$OMP CRITICAL(dos_)
        if (dosssum) then
          dt(:,1)=dt(:,1)-g(:)
        else
          dt(:,ispn)=dt(:,ispn)-g(:)
        end if
!$OMP END CRITICAL(dos_)
      end do
!$OMP END DO
      deallocate(f,g)
!$OMP END PARALLEL
      call freethd(nthd)
    end do
! output to file
    write(fname,'("PDOS_S",I2.2,"_A",I4.4)') is,ia
    open(50,file=trim(fname)//trim(fext),form='FORMATTED')
    do ispn=1,nsd
      do l=l0,l1
        do iw=1,nwplot
          write(50,'(2G18.10)') w(iw),dp(iw,l,ispn)*sps(ispn)
        end do
        write(50,'("     ")')
      end do
    end do
    close(50)
  end do
end do
!------------------------------------------!
!     irreducible representations file     !
!------------------------------------------!
if (lmirep) then
  open(50,file='ELMIREP'//trim(fext),form='FORMATTED')
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      do l=0,lmaxdos
        do m=-l,l
          lm=idxlm(l,m)
          write(50,'(" l = ",I2,", m = ",I2,", lm= ",I3," : ",G18.10)') l,m, &
           lm,elm(lm,ias)
        end do
      end do
    end do
  end do
  close(50)
end if
!--------------------------!
!     interstitial DOS     !
!--------------------------!
open(50,file='IDOS'//trim(fext),form='FORMATTED')
do ispn=1,nsd
  do iw=1,nwplot
    write(50,'(2G18.10)') w(iw),dt(iw,ispn)*sps(ispn)
  end do
  write(50,'("     ")')
end do
close(50)
! write the total DOS to test file
call writetest(10,'total DOS',nv=nwplot*nsd,tol=2.d-2,rva=dt)
deallocate(bc,sc,w,e,dt,dp)
if (lmirep) deallocate(elm,ulm)
return
end subroutine
!EOC

