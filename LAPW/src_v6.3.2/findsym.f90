
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsym
! !INTERFACE:
subroutine findsym(apl1,apl2,nsym,lspl,lspn,iea)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   apl1 : first set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   apl2 : second set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   nsym : number of symmetries (out,integer)
!   lspl : spatial rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   lspn : spin rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   iea  : equivalent atom index for each symmetry
!          (out,integer(iea(natmmax,nspecies,48))
! !DESCRIPTION:
!   Finds the symmetries which rotate one set of atomic positions into another.
!   Both sets of positions differ only by a translation vector and have the same
!   muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
!   symmetry element consists of a spatial rotation of the atomic position
!   vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
!   the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
!   returned as indices of elements in the Bravais lattice point group. An
!   index to equivalent atoms is stored in the array {\tt iea}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!   Fixed use of proper rotations for spin, February 2008 (L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: apl1(3,maxatoms,maxspecies)
real(8), intent(in) :: apl2(3,maxatoms,maxspecies)
integer, intent(out) :: nsym
integer, intent(out) :: lspl(48)
integer, intent(out) :: lspn(48)
integer, intent(out) :: iea(natmmax,nspecies,48)
! local variables
integer isym,jsym,jsym0,jsym1
integer is,ia,ias,ja,jas,md,n
real(8) sl(3,3),sc(3,3),v(3),t1
! automatic arrays
integer jea(natmmax,nspecies)
real(8) apl3(3,natmmax)
! allocatable arrays
complex(8), allocatable :: dmat(:,:,:,:)
! external functions
real(8) dnrm2
external dnrm2
nsym=0
! loop over lattice symmetries (spatial rotations)
do isym=1,nsymlat
! make real copy of lattice rotation symmetry
  sl(:,:)=dble(symlat(:,:,isym))
! loop over species
  do is=1,nspecies
! map apl1 coordinates to [0,1) and store in apl3
    do ia=1,natoms(is)
      apl3(:,ia)=apl1(:,ia,is)
      call r3frac(epslat,apl3(:,ia))
    end do
    do ja=1,natoms(is)
! apply lattice symmetry to atomic positions
      v(:)=sl(:,1)*apl2(1,ja,is)+sl(:,2)*apl2(2,ja,is)+sl(:,3)*apl2(3,ja,is)
! map coordinates to [0,1)
      call r3frac(epslat,v)
! check if atomic positions are invariant
      do ia=1,natoms(is)
        t1=abs(apl3(1,ia)-v(1))+abs(apl3(2,ia)-v(2))+abs(apl3(3,ia)-v(3))
        if (t1.lt.epslat) then
! equivalent atom index
          jea(ia,is)=ja
          goto 10
        end if
      end do
! not invariant so try new spatial rotation
      goto 40
10 continue
    end do
  end do
! all atomic positions invariant at this point
  jsym=1
! spin polarised case
  if (spinpol) then
! check invariance of magnetic fields under global spin rotation
    if (spinorb) then
! with spin-orbit coupling spin rotation equals spatial rotation
      jsym0=isym
      jsym1=isym
    else
! without spin-orbit coupling spin rotation independent of spatial rotation
      jsym0=1
      jsym1=nsymlat
    end if
    do jsym=jsym0,jsym1
! determinant of the symmetry matrix
      md=symlatd(jsym)
      sc(:,:)=dble(md)*symlatc(:,:,jsym)
! rotate global field and check invariance using proper part of symmetry matrix
      v(:)=sc(:,1)*bfieldc0(1)+sc(:,2)*bfieldc0(2)+sc(:,3)*bfieldc0(3)
      t1=abs(bfieldc0(1)-v(1))+abs(bfieldc0(2)-v(2))+abs(bfieldc0(3)-v(3))
! if not invariant try a different global spin rotation
      if (t1.gt.epslat) goto 20
! rotate muffin-tin magnetic fields and check invariance
      do is=1,nspecies
        do ia=1,natoms(is)
! equivalent atom
          ja=jea(ia,is)
          v(:)=sc(:,1)*bfcmt0(1,ja,is) &
              +sc(:,2)*bfcmt0(2,ja,is) &
              +sc(:,3)*bfcmt0(3,ja,is)
          t1=abs(bfcmt0(1,ia,is)-v(1)) &
            +abs(bfcmt0(2,ia,is)-v(2)) &
            +abs(bfcmt0(3,ia,is)-v(3))
! if not invariant try a different global spin rotation
          if (t1.gt.epslat) goto 20
        end do
      end do
! all fields invariant
      goto 30
20 continue
! end loop over global spin rotations
    end do
! magnetic fields not invariant so try different spatial rotation
    goto 40
  end if
30 continue
! check invariance of density matrices for fixed tensor moment calculations
  if (ftmtype.ne.0) then
    allocate(dmat(lmmaxdm,nspinor,lmmaxdm,nspinor))
    n=2*(lmmaxdm*nspinor)**2
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! equivalent atom
        ja=jea(ia,is)
        jas=idxas(ja,is)
! rotate the fixed tensor moment density matrix
        dmat(:,:,:,:)=0.d0
        call rotdmat(symlatc(:,:,isym),symlatc(:,:,jsym),lmaxdm,nspinor, &
         lmmaxdm,dmftm(:,:,:,:,jas),dmat)
! check invariance
        call daxpy(n,-1.d0,dmftm(:,:,:,:,ias),1,dmat,1)
        t1=dnrm2(n,dmat,1)/dble(n)
        if (t1.gt.epslat) then
          deallocate(dmat)
          goto 40
        end if
      end do
    end do
    deallocate(dmat)
  end if
! everything invariant so add symmetry to set
  nsym=nsym+1
  lspl(nsym)=isym
  lspn(nsym)=jsym
  do is=1,nspecies
    do ia=1,natoms(is)
      iea(ia,is,nsym)=jea(ia,is)
    end do
  end do
40 continue
! end loop over spatial rotations
end do
return
end subroutine
!EOC

