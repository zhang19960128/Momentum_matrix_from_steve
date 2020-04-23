MODULE io_espresso
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: read_es_header_qe, read_eigens, read_qekpt
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_es_header_qe()
  ! --------------------------------------------------------------------------
    !
    USE control, ONLY : kpts_type, dir => dirname
    USE es_tools, ONLY : es_header, lattice
    USE io_global, ONLY : stdout, iuinput
    USE iotk_module
    !
    IMPLICIT NONE
    !
    TYPE(es_header), POINTER :: read_es_header_qe
    !
    ! Cell parameters
    REAL(kind=8) :: a, at(3,3)
    ! Symmetry parameters
    INTEGER :: nsym
    LOGICAL :: linvsym, timeresym, t_resym
    INTEGER, ALLOCATABLE :: sym(:,:,:)
    REAL(kind=8), ALLOCATABLE :: latsym(:,:)
    ! Spin parameters
    LOGICAL :: lsda, sporb
    ! K-point parameters
    INTEGER :: nkpts, kptdim(3), kptoff(3)
    REAL(kind=8), ALLOCATABLE :: kpts_cart(:,:)
    ! Band structure params
    INTEGER :: nelec, nbnd, nspinor, fftdim(3)
    REAL(kind=8) :: nelec_temp
    !
    CHARACTER(iotk_attlenx) :: attr
    INTEGER :: i, ierr
    ! ------------------------------------------------------------------------
    !
    WRITE(stdout,'(A)') "Reading es_header ..."
    !
    CALL check_data_file(dir)
    !
! ----------------------------------------------------------------------------
!   ... Read the header data from file ...
! ----------------------------------------------------------------------------
    !
    ierr = 0
    CALL iotk_open_read(iuinput, FILE=trim(dir)//'.save'//'/'//"data-file.xml", IERR=ierr )
    IF ( ierr > 0 ) CALL call_error("read_es_header_qe", "Error opening xml file")
    !
    ! ... Cell parameters ...
    !
    CALL iotk_scan_begin( iuinput, "CELL" )
        CALL iotk_scan_dat(   iuinput, "LATTICE_PARAMETER", a )
        CALL iotk_scan_begin( iuinput, "DIRECT_LATTICE_VECTORS" )
            CALL iotk_scan_dat(   iuinput, "a1", at(:,1) )
            CALL iotk_scan_dat(   iuinput, "a2", at(:,2) )
            CALL iotk_scan_dat(   iuinput, "a3", at(:,3) )
        CALL iotk_scan_end(   iuinput, "DIRECT_LATTICE_VECTORS" )
    CALL iotk_scan_end( iuinput, "CELL" )
    !
    ! ... Symmetries ...
    !
    CALL iotk_scan_begin( iuinput, "SYMMETRIES" )
        CALL iotk_scan_dat( iuinput, "NUMBER_OF_SYMMETRIES", nsym )
        CALL iotk_scan_dat( iuinput, "INVERSION_SYMMETRY", linvsym )
        CALL iotk_scan_dat( iuinput, "DO_NOT_USE_TIME_REVERSAL", timeresym )
        CALL iotk_scan_dat( iuinput, "TIME_REVERSAL_FLAG", t_resym )
        !
        ALLOCATE( sym(3,3,nsym), latsym(3,nsym) )
        !
        DO i=1,nsym
            CALL iotk_scan_begin( iuinput, "SYMM"// TRIM( iotk_index( i ) ) )
                CALL iotk_scan_dat( iuinput, "ROTATION", sym(:,:,i))
                   sym(:,:,i)=transpose(sym(:,:,i))
                CALL iotk_scan_dat( iuinput, "FRACTIONAL_TRANSLATION", latsym(:,i))
            CALL iotk_scan_end( iuinput, "SYMM"// TRIM( iotk_index( i ) ) )
        END DO
    CALL iotk_scan_end( iuinput, "SYMMETRIES" )
    !
    ! ... Spin ...
    !
    CALL iotk_scan_begin (iuinput, "SPIN" )
        CALL iotk_scan_dat (iuinput, "LSDA", lsda)
        CALL iotk_scan_dat (iuinput, "SPIN-ORBIT_CALCULATION", sporb)
    CALL iotk_scan_end (iuinput, "SPIN" )
    !
    ! ... K-point data ...
    !
    CALL iotk_scan_begin( iuinput, "BRILLOUIN_ZONE" )
        CALL iotk_scan_dat( iuinput, "NUMBER_OF_K-POINTS", nkpts)
        !
        CALL iotk_scan_empty( iuinput, "MONKHORST_PACK_GRID", attr )
            CALL iotk_scan_attr(attr, "nk1",kptdim(1))
            CALL iotk_scan_attr(attr, "nk2",kptdim(2))
            CALL iotk_scan_attr(attr, "nk3",kptdim(3))
            !
        CALL iotk_scan_empty( iuinput, "MONKHORST_PACK_OFFSET", attr )
            CALL iotk_scan_attr(attr, "k1",kptoff(1))
            CALL iotk_scan_attr(attr, "k2",kptoff(2))
            CALL iotk_scan_attr(attr, "k3",kptoff(3))
            !
        ALLOCATE( kpts_cart(3,nkpts) )
        DO i=1,nkpts
            CALL iotk_scan_empty( iuinput, "K-POINT"//TRIM(iotk_index(i)), attr )
                CALL iotk_scan_attr(attr, "XYZ",kpts_cart(:,i))
        END DO
    CALL iotk_scan_end( iuinput, "BRILLOUIN_ZONE" )
    !
    ! non monkhorst pack grid currently implemented only for nxnxn
    IF ( trim(kpts_type) .eq. 'nonmp' ) THEN
      kptdim(1) = nint((1.0*nkpts)**(1.0/3.0))
      kptdim(2) = nint((1.0*nkpts)**(1.0/3.0))
      kptdim(3) = nint((1.0*nkpts)**(1.0/3.0))
    END IF
    !
    IF ( trim(kpts_type) .ne. 'path' ) THEN
      IF (any(kptdim == 0)) then
        CALL call_warning("KPT GRID is zero, ","go to data-file.xml and change BRILLOUIN_ZONE//MONKHORST_PACK_GRID")
        CALL call_error("read_es_header_qe","Please indicate kpt grid in data-file.xml")
      END IF
    END IF
    !
    ! ... Band structure ...
    !
    CALL iotk_scan_begin( iuinput, "BAND_STRUCTURE_INFO" )
        CALL iotk_scan_dat( iuinput, "NUMBER_OF_SPIN_COMPONENTS", nspinor)
        CALL iotk_scan_dat( iuinput, "NUMBER_OF_BANDS", nbnd)
        CALL iotk_scan_dat( iuinput, "NUMBER_OF_ELECTRONS", nelec_temp)
    CALL iotk_scan_end( iuinput, "BAND_STRUCTURE_INFO" )
    nelec = int(nelec_temp)
    !
    CALL iotk_scan_begin( iuinput, "PLANE_WAVES" )
        CALL iotk_scan_empty(iuinput, "FFT_GRID", attr)
            CALL iotk_scan_attr(attr, "nr1",fftdim(1))
            CALL iotk_scan_attr(attr, "nr2",fftdim(2))
            CALL iotk_scan_attr(attr, "nr3",fftdim(3))
    CALL iotk_scan_end( iuinput, "PLANE_WAVES")
    !
    CALL iotk_close_read( iuinput )
    !
! ----------------------------------------------------------------------------
!   ... Save the header to data structure ...
! ----------------------------------------------------------------------------
    !
    ALLOCATE( read_es_header_qe )
    !
    ! ... Lattice and lattice symmetries ...
    read_es_header_qe%lattice = lattice(at, sym, latsym)
    !
    ! ... Other Symmetries ...
    IF ( trim(kpts_type).eq.'nonmp' .or. trim(kpts_type).eq.'path' ) THEN
      read_es_header_qe%gen_time_reversal = .false.
    ELSE
      read_es_header_qe%gen_time_reversal = (.not. timeresym) .and. t_resym
    END IF
    !
    ! ... Spin ...
    read_es_header_qe%spinorb = sporb
    read_es_header_qe%nspinor = nspinor
    IF (lsda) THEN
      read_es_header_qe%nspin = 2
    ELSE
      read_es_header_qe%nspin = 1
    END IF
    !
    ! ... K-point data ...
    read_es_header_qe%nkpt = nkpts
    read_es_header_qe%kptdim = kptdim
    !
    forall(i=1:3) read_es_header_qe%kptoffset(i) = kptoff(i)/(2d0*kptdim(i))
    !
    ! kcoord use conventional fractional k space coordinates
    ! kptns is the crystal coord
    ALLOCATE( read_es_header_qe%kptns(3,nsym,nkpts) )
    read_es_header_qe%kptns(:,1,:) = matmul(transpose(at),kpts_cart/a)
    !
    ! ... Electronic info ...
    read_es_header_qe%nelec = nelec
    read_es_header_qe%nband = nbnd
    read_es_header_qe%fftdim = fftdim
    !
    ! ... Echo back ...
    WRITE(stdout,*)
    WRITE(stdout,'(A, L  )') 'LSDA            = ', lsda
    WRITE(stdout,'(A,3I3 )') 'K-grid          = ', read_es_header_qe%kptdim
    WRITE(stdout,'(A, I8 )') 'Num of kpoints  = ', read_es_header_qe%nkpt
    WRITE(stdout,'(A, I3 )') 'Num of sym      = ', read_es_header_qe%lattice%nsym
    WRITE(stdout,'(A, I5 )') 'Num of elecs    = ', read_es_header_qe%nelec
    WRITE(stdout,'(A,3I5 )') 'FFT-grid        = ', read_es_header_qe%fftdim
    WRITE(stdout,'(A, L  )') 'Time Re. Sym.   = ', read_es_header_qe%gen_time_reversal
    !
  END FUNCTION read_es_header_qe
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE check_data_file(dir)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: dir
    LOGICAL :: datafile
    !
    INQUIRE(FILE=trim(dir)//'.save'//'/'//"data-file.xml", EXIST=datafile)
    !
    IF (.not. datafile) CALL call_error('check_data_file', "data-file.xml is missing")
  END SUBROUTINE
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_eigens( prefix, header_kl, output_type )
  ! --------------------------------------------------------------------------
    !
    USE constants, ONLY : dp, hartree_ev
    USE es_tools, ONLY : es_header
    USE io_global, ONLY : stdout, iuinput
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(inout) :: header_kl
    CHARACTER(LEN=*), INTENT(in) :: prefix
    INTEGER, INTENT(in) :: output_type
    !
    INTEGER :: nelec, nspin, nband, nkpt
    INTEGER :: ikpt, iband, ispin
    INTEGER :: band_ind, k_ind
    INTEGER :: occband, occ_multi
    REAL(dp) :: occ
    CHARACTER(LEN=1) :: is
    !
    IF (output_type==1) WRITE(stdout,'(A)') "QE fillings are determined by input!"
    IF (output_type==2) WRITE(stdout,'(A)') "Wannier fillings are determined by number of electrons! Only for insulators!"
    !
    nkpt = header_kl%nkpt
    nelec = header_kl%nelec
    nspin = header_kl%nspin
    nband = header_kl%nband
    !
    ALLOCATE ( header_kl%energy(nband,nspin,nkpt), header_kl%filling(nband,nspin,nkpt) )
    !
    DO ispin=1,nspin
      WRITE(is, fmt='(I1)') ispin
      OPEN(iuinput, file=trim(prefix) // '.unk' // '/' // trim(prefix)//'.eigen.'//is, status='old')
      !
      DO ikpt=1,nkpt
        DO iband=1,nband
          ! QE input, read fillings
          IF (output_type==1) read(iuinput,*) band_ind, k_ind, header_kl%energy(iband,ispin,ikpt), header_kl%filling(iband,ispin,ikpt)
          ! W90 input, decide fillings later
          IF (output_type==2) read(iuinput,*) band_ind, k_ind, header_kl%energy(iband,ispin,ikpt)
        END DO
      END DO
      !
    END DO
    !
    ! Decide occupation correction
    occ_multi = 3-nspin
    IF (header_kl%spinorb) occ_multi = 1
    occ = dble(occ_multi)
    !
    IF (output_type==1) THEN
      ! QE input, correct fillings for non-spin polarized calc
      header_kl%filling = header_kl%filling*occ
    END IF
    !
    IF (output_type==2) THEN
      ! written energy for QE has unit of Hartree, but wannier has unit of eV
      header_kl%energy = header_kl%energy/hartree_ev
      !
      ! W90 input, decide fillings based on electron count
      occband=int(nelec/occ)
      IF (mod(int(nelec),2)/=0) CALL call_error("read_eigens", "Num of elec is not even; ferromagnetism not implented!")
      header_kl%filling=0.0
      forall(iband=1:occband,ispin=1:nspin,ikpt=1:nkpt) header_kl%filling(iband,ispin,ikpt) = occ
    END IF
    !
  END SUBROUTINE read_eigens
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_qekpt( ikpt, kcoord, energy, filling, nband, ispin, tr, sporb )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp, dpc
    USE control, ONLY : dir => dirname
    USE es_tools, ONLY : kpt_data, is_tr, kpt_init, kpt_init_tr
    USE io_w90, ONLY : unk_kpoint_file, unk_gvect_file
    !
    IMPLICIT NONE
    !
    TYPE(kpt_data) :: read_qekpt
    INTEGER, INTENT(in) :: ikpt, ispin, nband
    REAL(dp), INTENT(in) :: kcoord(3)
    REAL(dp), INTENT(in) :: energy(:), filling(:)
    LOGICAL, INTENT(in) :: tr, sporb
    !
    CHARACTER(LEN=256) :: kdir, kdir2, gdir
    INTEGER, ALLOCATABLE :: g(:,:)
    COMPLEX(dpc), ALLOCATABLE :: bands(:,:,:)
    !
    IF (sporb) THEN
      kdir = unk_kpoint_file(trim(dir) // '.unk', ikpt, ispin)
      kdir2 = unk_kpoint_file(trim(dir) // '.unk', ikpt, 2)
      CALL read_bands(kdir, bands, kdir2)
      !
    ELSE
      kdir = unk_kpoint_file(trim(dir) // '.unk', ikpt, ispin)
      CALL read_bands(kdir, bands)
      !
    END IF
    !
    gdir = unk_gvect_file(trim(dir) // '.unk', ikpt, 1)
    CALL read_gvect(gdir, g)
    !
    IF (is_tr(kcoord)) THEN
      CALL kpt_init_tr(kcoord, read_qekpt, g, bands, energy, filling, tr)
    ELSE
      CALL kpt_init(kcoord, read_qekpt, g, bands, energy, filling)
    END IF
    !
    IF (allocated(g)) deallocate(g)
    IF (allocated(bands)) deallocate(bands)
    !
  END FUNCTION read_qekpt
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_bands( filename, bands, filename2 )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dpc
    USE io_global, ONLY : iuinput
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: filename2
    COMPLEX(dpc), ALLOCATABLE, INTENT(out) :: bands(:,:,:)
    !
    INTEGER :: npw, nspinor, nband, npw2, nband2
    INTEGER :: ipw, ikpt, iband, ikpt2
    !
    nspinor = 1
    IF (present(filename2)) nspinor = 2
    !
    OPEN( iuinput, FILE = trim(filename), FORM='unformatted' )
    READ( iuinput ) npw, ikpt, nband
    !
    ALLOCATE ( bands( npw, nspinor, nband ) )
    !
    DO iband=1,nband
      READ( iuinput ) bands(:,1,iband)
    END DO
    CLOSE( iuinput )
    !
    IF (present(filename2)) THEN
      !
      OPEN( iuinput, FILE = trim(filename2), FORM='unformatted' )
      READ( iuinput ) npw2, ikpt2, nband2
      !
      IF (npw2/=npw .or. ikpt2/=ikpt .or. nband2/=nband) CALL call_error("read_bands", 'spin orbital calculation wfc not consistent!')
      !
      DO iband=1,nband
        READ( iuinput ) bands(:,2,iband)
      END DO
      CLOSE( iuinput )
      !
    END IF
    !
  END SUBROUTINE read_bands
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_gvect ( filename, g )
  ! --------------------------------------------------------------------------
    USE io_global, ONLY : iuinput
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER, ALLOCATABLE, INTENT(out) :: g(:,:)
    !
    INTEGER :: ipw, npw
    !
    OPEN( iuinput, FILE = TRIM(filename), FORM='unformatted' )
    READ( iuinput ) npw
    !
    ALLOCATE ( g( 3, npw ) )
    !
    DO ipw=1,npw
      READ( iuinput ) g(:,ipw)
    END DO
    !
    CLOSE(iuinput)
    !
  END SUBROUTINE read_gvect
  !
END MODULE io_espresso
