SUBROUTINE read_input()
  USE modmain
  USE modrandom, ONLY: randomu
  IMPLICIT NONE 
  ! local variables
  logical :: lv
  INTEGER :: is,ia,ias,ios
  INTEGER :: i,j,k,l,n,p
  REAL(8) :: sc,sc1,sc2,sc3
  REAL(8) :: scu,scu1,scu2,scu3
  REAL(8) :: solscf,zn,a,b
  REAL(8) :: axang(4),rot(3,3)
  REAL(8) :: v1(3),v2(3),t1
  character(256) :: block,symb,str
    
  CALL default_apwlo()
  CALL default_atoms()
  CALL default_charge_moment_current()
  CALL default_convergence()
  CALL default_core_states()

  CALL default_density_pot_xc()
  CALL default_dos_optics_response()
  CALL default_electric_vector_pot()
  CALL default_force_stress()
  CALL default_gkvectors()
  
  CALL default_gvectors()
  CALL default_hamiltonian()
  CALL default_kpoints()
  CALL default_lattice()
  CALL default_misc()
  
  CALL default_mixing()
  CALL default_muffin_tins()
  CALL default_oep_hf()
  CALL default_plotting()
  CALL default_qpoints()
  
  CALL default_spin()
  CALL default_states()
  CALL default_symmetry()

  sc=1.d0
  sc1=1.d0
  sc2=1.d0
  sc3=1.d0
  solscf=1.d0
  axang(:)=0.d0

  !--------------------------!
  !     read from elk.in     !
  !--------------------------!
  open(50,file='elk.in',status='OLD',form='FORMATTED',iostat=ios)
  IF(ios.ne.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): error opening elk.in")')
    WRITE(*,*)
    STOP 
  ENDIF 
  10 continue
  
  read(50,*,end=30) block
  ! check for a comment
  IF((scan(trim(block),'!')==1) .or. (scan(trim(block),'#')==1)) GOTO 10
  SELECT CASE(trim(block))
  CASE('tasks')
    DO i=1,maxtasks
      read(50,'(A256)',err=20) str
      IF(trim(str).eq.'') THEN 
        IF(i.eq.1) THEN 
          WRITE(*,*)
          WRITE(*,'("Error(readinput): no tasks to perform")')
          WRITE(*,*)
          STOP 
        ENDIF 
        ntasks=i-1
        GOTO 10
      ENDIF 
      read(str,*,iostat=ios) tasks(i)
      IF(ios.ne.0) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(readinput): error reading tasks")')
        WRITE(*,'("(blank line required after tasks block)")')
        WRITE(*,*)
        STOP 
      ENDIF 
    ENDDO 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): too many tasks")')
    WRITE(*,'("Adjust maxtasks in modmain and recompile code")')
    WRITE(*,*)
    STOP 
  !CASE('species')
  !  ! generate a species file
  !  CALL genspecies(50)
  !CASE('fspecies')
  !  ! generate fractional species files
  !  DO is=1,maxspecies
  !    read(50,'(A256)',err=20) str
  !    IF(trim(str).eq.'') GOTO 10
  !    read(str,*,iostat=ios) zn,symb
  !    IF(zn.gt.-1.d0+epsocc) THEN 
  !      WRITE(*,*)
  !      WRITE(*,'("Error(readinput): fractional nuclear Z > -1 : ",G18.10)') zn
  !      WRITE(*,*)
  !      STOP 
  !    ENDIF 
  !    CALL genfspecies(zn,symb)
  !  ENDDO 
  !  WRITE(*,*)
  !  WRITE(*,'("Error(readinput): too many fractional nucleus species")')
  !  WRITE(*,*)
  !  STOP 
  CASE('avec')
    read(50,*,err=20) avec(:,1)
    read(50,*,err=20) avec(:,2)
    read(50,*,err=20) avec(:,3)
  CASE('primcell')
    read(50,*,err=20) primcell
  CASE('tshift')
    read(50,*,err=20) tshift
  CASE('autokpt')
    read(50,*,err=20) autokpt
  CASE('radkpt')
    read(50,*,err=20) radkpt
    IF(radkpt.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): radkpt <= 0 : ",G18.10)') radkpt
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('ngridk')
    read(50,*,err=20) ngridk(:)
    IF((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
      WRITE(*,*)
      STOP 
    ENDIF 
    autokpt=.false.
  CASE('vkloff')
    read(50,*,err=20) vkloff(:)
  CASE('reducek')
    read(50,*,err=20) reducek
  CASE('ngridq')
    read(50,*,err=20) ngridq(:)
    IF((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('reduceq')
    read(50,*,err=20) reduceq
  CASE('rgkmax')
    read(50,*,err=20) rgkmax
    IF(rgkmax.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('gmaxvr')
    read(50,*,err=20) gmaxvr
  CASE('lmaxapw')
    read(50,*,err=20) lmaxapw
    IF(lmaxapw.lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): lmaxapw < 0 : ",I8)') lmaxapw
      WRITE(*,*)
      STOP 
    ENDIF 
    IF(lmaxapw.ge.maxlapw) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): lmaxapw too large : ",I8)') lmaxapw
      WRITE(*,'("Adjust maxlapw in modmain and recompile code")')
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('lmaxo','lmaxvr')
    read(50,*,err=20) lmaxo
    IF(lmaxo.lt.3) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): lmaxo < 3 : ",I8)') lmaxo
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('lmaxi','lmaxinr')
    read(50,*,err=20) lmaxi
    IF(lmaxi.lt.1) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): lmaxi < 1 : ",I8)') lmaxi
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('lmaxmat')
    read(50,*,err=20)
    WRITE(*,*)
    WRITE(*,'("Info(readinput): variable ''lmaxmat'' is no longer used")')
  CASE('fracinr')
    read(50,*,err=20) fracinr
  CASE('trhonorm')
    read(50,*,err=20) trhonorm
  CASE('spinpol')
    read(50,*,err=20) spinpol
  CASE('spinorb')
    read(50,*,err=20) spinorb
  CASE('socscf')
    read(50,*,err=20) socscf
    IF(socscf.lt.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): socscf < 0 : ",G18.10)') socscf
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('xctype')
    read(50,'(A256)',err=20) str
    str=trim(str)//' 0 0'
    read(str,*,err=20) xctype
  CASE('xctsp')
    read(50,'(A256)',err=20) str
    str=trim(str)//' 0 0'
    read(str,*,err=20) xctsp
  CASE('stype')
    read(50,*,err=20) stype
  CASE('swidth')
    read(50,*,err=20) swidth
    IF(swidth.lt.1.d-9) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
       swidth
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('autoswidth')
    read(50,*,err=20) autoswidth
  CASE('mstar')
    read(50,*,err=20) mstar
    IF(mstar.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): mstar <= 0 : ",G18.10)') mstar
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('epsocc')
    read(50,*,err=20) epsocc
    IF(epsocc.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('epschg')
    read(50,*,err=20) epschg
    IF(epschg.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('nempty','nempty0')
    read(50,*,err=20) nempty0
    IF(nempty0.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): nempty <= 0 : ",G18.10)') nempty0
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('mixtype')
    read(50,*,err=20) mixtype
  CASE('amixpm','beta0','betamax')
    IF(trim(block).eq.'amixpm') THEN 
      read(50,*,err=20) amixpm(:)
    ELSEIF(trim(block).eq.'beta0') THEN 
      read(50,*,err=20) amixpm(1)
    else
      read(50,*,err=20) amixpm(2)
    ENDIF 
    IF(amixpm(1).lt.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): beta0 [amixpm(1)] < 0 : ",G18.10)') amixpm(1)
      WRITE(*,*)
      STOP 
    ENDIF 
    IF((amixpm(2).lt.0.d0).or.(amixpm(2).gt.1.d0)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): betamax [amixpm(2)] not in [0,1] : ",G18.10)')&
       amixpm(2)
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('mixsdb')
    read(50,*,err=20) mixsdb
    IF(mixsdb.lt.2) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): mixsdb < 2 : ",I8)') mixsdb
      WRITE(*,*)
      STOP 
    ENDIF 
CASE('broydpm')
  read(50,*,err=20) broydpm(:)
  IF((broydpm(1).lt.0.d0).or.(broydpm(1).gt.1.d0).or. &
      (broydpm(2).lt.0.d0).or.(broydpm(2).gt.1.d0)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): invalid Broyden mixing parameters : ",&
     &2G18.10)') broydpm
    WRITE(*,*)
    STOP 
  ENDIF 
CASE('maxscl')
  read(50,*,err=20) maxscl
  IF(maxscl.lt.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    WRITE(*,*)
    STOP 
  ENDIF 
CASE('epspot')
  read(50,*,err=20) epspot
CASE('epsengy')
  read(50,*,err=20) epsengy
CASE('epsforce')
  read(50,*,err=20) epsforce
CASE('epsstress')
  read(50,*,err=20) epsstress
CASE('sppath')
  read(50,*,err=20) sppath
  sppath=adjustl(sppath)
CASE('scrpath')
  read(50,*,err=20) scrpath
CASE('molecule')
  read(50,*,err=20) molecule
CASE('atoms')
  read(50,*,err=20) nspecies
  IF(nspecies.le.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
    WRITE(*,*)
    STOP 
  ENDIF 
  IF(nspecies.gt.maxspecies) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
    WRITE(*,'("Adjust maxspecies in modmain and recompile code")')
    WRITE(*,*)
    STOP 
  ENDIF 
  DO is=1,nspecies
    read(50,*,err=20) spfname(is)
    spfname(is)=adjustl(spfname(is))
    read(50,*,err=20) natoms(is)
    IF(natoms(is).le.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,*)
      STOP 
    ENDIF 
    IF(natoms(is).gt.maxatoms) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): natoms too large : ",I8)') natoms(is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'("Adjust maxatoms in modmain and recompile code")')
      WRITE(*,*)
      STOP 
    ENDIF 
    DO ia=1,natoms(is)
      read(50,'(A256)',err=20) str
      str=trim(str)//' 0.0 0.0 0.0'
      read(str,*,err=20) atposl(:,ia,is),bfcmt0(:,ia,is)
    ENDDO 
  ENDDO 
  CASE('dosocc')
    read(50,*,err=20) dosocc
  CASE('dosmsum')
    read(50,*,err=20) dosmsum
  CASE('dosssum')
    read(50,*,err=20) dosssum
  CASE('lmirep')
    read(50,*,err=20) lmirep
  CASE('atpopt')
    read(50,*,err=20) atpopt
  CASE('maxatpstp','maxatmstp')
    read(50,*,err=20) maxatpstp
    IF(maxatpstp.le.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): maxatpstp <= 0 : ",I8)') maxatpstp
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('tau0atp','tau0atm')
    read(50,*,err=20) tau0atp
  CASE('deltast')
    read(50,*,err=20) deltast
    IF(deltast.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): deltast <= 0 : ",G18.10)') deltast
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('latvopt')
    read(50,*,err=20) latvopt
  CASE('maxlatvstp')
    read(50,*,err=20) maxlatvstp
    IF(maxlatvstp.le.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): maxlatvstp <= 0 : ",I8)') maxlatvstp
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('tau0latv')
    read(50,*,err=20) tau0latv
  CASE('nstfsp')
    read(50,*,err=20)
    WRITE(*,*)
    WRITE(*,'("Info(readinput): variable ''nstfsp'' is no longer used")')
  CASE('lradstp')
    read(50,*,err=20) lradstp
    IF(lradstp.le.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('chgexs')
    read(50,*,err=20) chgexs
  CASE('nprad')
    read(50,*,err=20)
    WRITE(*,*)
    WRITE(*,'("Info(readinput): variable ''nprad'' is no longer used")')
  CASE('scissor')
    read(50,*,err=20) scissor
CASE('optcomp')
  DO i=1,27
    read(50,'(A256)',err=20) str
    IF(trim(str).eq.'') THEN 
      IF(i.eq.1) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(readinput): empty optical component list")')
        WRITE(*,*)
        STOP 
      ENDIF 
      noptcomp=i-1
      GOTO 10
    ENDIF 
    str=trim(str)//' 1 1'
    read(str,*,iostat=ios) optcomp(:,i)
    IF(ios /= 0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): error reading optical component list")')
      WRITE(*,'("(blank line required after optcomp block)")')
      WRITE(*,*)
      STOP 
    ENDIF 
    IF((optcomp(1,i).lt.1).or.(optcomp(1,i).gt.3).or. &
        (optcomp(2,i).lt.1).or.(optcomp(2,i).gt.3).or. &
        (optcomp(3,i).lt.1).or.(optcomp(3,i).gt.3)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): invalid optcomp : ",3I8)') optcomp
      WRITE(*,*)
      STOP 
    ENDIF 
  ENDDO 
  WRITE(*,*)
  WRITE(*,'("Error(readinput): optical component list too long")')
  WRITE(*,*)
  STOP 
  CASE('autormt')
    read(50,*,err=20)
    WRITE(*,*)
    WRITE(*,'("Info(readinput): variable ''autormt'' is no longer used")')
  CASE('rmtdelta')
    read(50,*,err=20) rmtdelta
    IF(rmtdelta.lt.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Warning(readinput): rmtdelta < 0 : ",G18.10)') rmtdelta
    ENDIF 
  CASE('isgkmax')
    read(50,*,err=20) isgkmax
  CASE('nosym')
    read(50,*,err=20) lv
    IF(lv) symtype=0
  CASE('symtype')
    read(50,*,err=20) symtype
    IF((symtype.lt.0).or.(symtype.gt.2)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): symtype not defined : ",I8)') symtype
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('notes')
    IF(allocated(notes)) DEALLOCATE(notes)
    ALLOCATE(notes(0))
    notelns=0
    DO 
      read(50,'(A80)') str
      IF(trim(str)=='') GOTO 10
      notelns=notelns+1
      CALL addstr(notes)
    ENDDO 
  CASE('tforce')
    read(50,*,err=20) tforce
  CASE('maxitoep')
    read(50,*,err=20) maxitoep
    IF(maxitoep.lt.1) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('kstlist')
    DO i=1,maxkst
      read(50,'(A256)',err=20) str
      IF(trim(str).eq.'') THEN 
        IF(i.eq.1) THEN 
          WRITE(*,*)
          WRITE(*,'("Error(readinput): empty k-point and state list")')
          WRITE(*,*)
          STOP 
        ENDIF 
        nkstlist=i-1
        GOTO 10
      ENDIF 
      read(str,*,iostat=ios) kstlist(:,i)
      IF(ios.ne.0) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(readinput): error reading k-point and state list")')
        WRITE(*,'("(blank line required after kstlist block)")')
        WRITE(*,*)
        STOP 
      ENDIF 
    ENDDO 
    WRITE(*,*)
    WRITE(*,'("Error(readinput): k-point and state list too long")')
    WRITE(*,*)
    STOP 
  CASE('vklem')
    read(50,*,err=20) vklem
  CASE('deltaem')
    read(50,*,err=20) deltaem
  CASE('ndspem')
    read(50,*,err=20) ndspem
    IF((ndspem.lt.1).or.(ndspem.gt.4)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('ptnucl')
    read(50,*,err=20) ptnucl
  CASE('tefvr','tseqr')
    read(50,*,err=20) tefvr
  CASE('tefvit','tseqit')
    read(50,*,err=20) tefvit
  CASE('minitefv','minseqit')
    read(50,*,err=20) minitefv
    IF(minitefv.lt.1) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): minitefv < 1 : ",I8)') minitefv
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('emaxelnes')
    read(50,*,err=20) emaxelnes
  CASE('wsfac')
    read(50,*,err=20) wsfac(:)
  CASE('hybrid')
    read(50,*,err=20) hybrid
  CASE('hybridc','hybmix')
    read(50,*,err=20) hybridc
    IF((hybridc.lt.0.d0).or.(hybridc.gt.1.d0)) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): invalid hybridc : ",G18.10)') hybridc
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('nxoapwlo','nxapwlo')
    read(50,*,err=20) nxoapwlo
    IF(nxoapwlo.lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): nxoapwlo < 0 : ",I8)') nxoapwlo
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('nxlo')
    read(50,*,err=20) nxlo
    IF(nxlo.lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): nxlo < 0 : ",I8)') nxlo
      WRITE(*,*)
      STOP 
    ENDIF 
  CASE('tempk')
    read(50,*,err=20) tempk
    IF(tempk.le.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readinput): tempk <= 0 : ",G18.10)') tempk
      WRITE(*,*)
      STOP 
    ENDIF 
    ! set Fermi-Dirac smearing
    stype=3
    ! set the smearing width
    swidth=kboltz*tempk
  CASE('rmtall')
    read(50,*,err=20) rmtall
  CASE('stable')
    read(50,*,err=20) lv
    IF(lv) THEN 
      nxoapwlo=max(nxoapwlo,1)
      mixtype=3
      broydpm(1)=min(broydpm(1),0.01d0)
      broydpm(2)=min(broydpm(2),0.04d0)
      msmooth=max(msmooth,8)
      WRITE(*,*)
      WRITE(*,'("Info(readinput): parameters set by stable option")')
      WRITE(*,'(" nxoapwlo : ",I4)') nxoapwlo
      WRITE(*,'(" mixtype : ",I4)') mixtype
      WRITE(*,'(" broydpm : ",2G18.10)') broydpm
      WRITE(*,'(" msmooth : ",I4)') msmooth
    ENDIF 
  !
  CASE('metagga')
    read(50,*,err=20) lv
    IF(lv) THEN 
      nempty0=max(nempty0,10.d0)
      lradstp=1
      nrmtscf=max(nrmtscf,2.d0)
      msmooth=max(msmooth,4)
      WRITE(*,*)
      WRITE(*,'("Info(readinput): parameters set by metagga option")')
      WRITE(*,'(" nempty0 : ",G18.10)') nempty0
      WRITE(*,'(" lradstp : ",I4)') lradstp
      WRITE(*,'(" nrmtscf : ",G18.10)') nrmtscf
      WRITE(*,'(" msmooth : ",I4)') msmooth
    ENDIF
  !
  CASE('autolinengy')
    read(50,*,err=20) autolinengy
  !
  CASE('')
    GOTO 10
  case default
    WRITE(*,*)
    WRITE(*,'("Error(readinput): invalid block name : ",A)') trim(block)
    WRITE(*,*)
    STOP 
  end select
  GOTO 10
  20 continue
  WRITE(*,*)
  WRITE(*,'("Error(readinput): error reading from elk.in")')
  WRITE(*,'("Problem occurred in ''",A,"'' block")') trim(block)
  WRITE(*,'("Check input convention in manual")')
  WRITE(*,*)
  STOP 
  30 continue
  close(50)

  ! scale the speed of light
  solsc=sol*solscf

  ! scale and rotate the lattice vectors (not referenced again in code)
  avec(:,1)=sc1*avec(:,1)
  avec(:,2)=sc2*avec(:,2)
  avec(:,3)=sc3*avec(:,3)
  avec(:,:)=sc*avec(:,:)
  t1=axang(4)
  IF(t1.ne.0.d0) THEN 
    t1=t1*pi/180.d0
    CALL axangrot(axang(:),t1,rot)
    DO i=1,3
      v1(:)=avec(:,i)
      CALL r3mv(rot,v1,avec(:,i))
    ENDDO 
  ENDIF 

  ! randomise lattice vectors if required
  IF(rndavec.gt.0.d0) THEN 
    DO i=1,3
      DO j=1,3
        t1=rndavec*(randomu()-0.5d0)
        avec(i,j)=avec(i,j)+t1
      ENDDO 
    ENDDO 
  ENDIF 

  ! case of isolated molecule
  IF(molecule) THEN 
    ! convert atomic positions from Cartesian to lattice coordinates
    CALL r3minv(avec,ainv)
    DO is=1,nspecies
      DO ia=1,natoms(is)
        CALL r3mv(ainv,atposl(:,ia,is),v1)
        atposl(:,ia,is)=v1(:)
      ENDDO 
    ENDDO 
  ENDIF 

  ! randomise atomic positions if required
  IF(rndatposc.gt.0.d0) THEN 
    CALL r3minv(avec,ainv)
    DO is=1,nspecies
      DO ia=1,natoms(is)
        CALL r3mv(avec,atposl(:,ia,is),v1)
        DO i=1,3
          t1=rndatposc*(randomu()-0.5d0)
          v1(i)=v1(i)+t1
        ENDDO 
        CALL r3mv(ainv,v1,atposl(:,ia,is))
      ENDDO 
    ENDDO 
  ENDIF 

  ! randomise the muffin-tin magnetic fields if required
  !IF(rndbfcmt > 0.d0) THEN 
  !  DO is=1,nspecies
  !    DO ia=1,natoms(is)
  !      DO i=1,3
  !        t1=rndbfcmt*(randomu()-0.5d0)
  !        bfcmt0(i,ia,is)=bfcmt0(i,ia,is)+t1
  !      ENDDO 
  !    ENDDO 
  !  ENDDO 
  !ENDIF 
  
  ! find primitive cell if required
  !IF(primcell) CALL findprimcell
  
  ! scale the ultracell vectors if required
  !avecu(:,1)=scu1*avecu(:,1)
  !avecu(:,2)=scu2*avecu(:,2)
  !avecu(:,3)=scu3*avecu(:,3)
  !avecu(:,:)=scu*avecu(:,:)

  ! read in atomic species data
  CALL read_species_files()
  
  WRITE(*,*) 'Finished reading input file'
  RETURN 

CONTAINS 

  SUBROUTINE addstr(slist)
    IMPLICIT NONE 
    ! arguments
    character(256), intent(inout), ALLOCATABLE :: slist(:)
    ! ALLOCATABLE arrays
    character(256), ALLOCATABLE :: stmp(:)
    n=size(slist)
    ALLOCATE(stmp(n))
    stmp(1:n)=slist(1:n)
    DEALLOCATE(slist)
    ALLOCATE(slist(n+1))
    slist(1:n)=stmp(1:n)
    slist(n+1)=str
    DEALLOCATE(stmp)
    RETURN 
  END SUBROUTINE 

END SUBROUTINE 