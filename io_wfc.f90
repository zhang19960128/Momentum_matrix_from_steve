MODULE io_wfc
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: trans_qewfc
CONTAINS
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE trans_qewfc(header, output_type)
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp, dpc
    USE control, ONLY : dir => dirname
    USE es_tools, ONLY : es_header
    USE io_global, ONLY : ionode, ionode_id, stdout, iuout
    USE io_w90, ONLY : unk_kpoint_file, unk_gvect_file
    USE mp_global, ONLY : nproc, proc_id, mp_barrier, mp_gather, mp_comm_all
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(in) :: header
    INTEGER, INTENT(in) :: output_type
    !
    INTEGER :: nkpt, nspinor, nbnd, nelec, nspin, ncolum
    INTEGER :: ikptset, kpt_ind
    INTEGER, ALLOCATABLE:: g(:,:), link_kpt(:,:)
    REAL(dp), ALLOCATABLE :: energy(:,:), filling(:,:)
    REAL(dp), ALLOCATABLE :: energy_all(:,:,:), filling_all(:,:,:)
    COMPLEX(dpc), ALLOCATABLE :: bands(:,:,:)
    CHARACTER(LEN=256) :: kdir, kpt_file, gvec_file, eig_file
    CHARACTER(LEN=1) :: is
        integer                                                              :: iband, ikpt,ispin, ipw

 
    IF (ionode) THEN
      CALL system('mkdir '//trim(dir)//'.unk')
      !
      WRITE(stdout,'(A)') '============================'
      WRITE(stdout,'(A)') '  Rewrite QE wavefunctions  '
      WRITE(stdout,'(A)') '  VERSION 1.1.0             '
    END IF
    !
    ! Get dimensions for calculation
    nkpt = header%nkpt
    nspinor = header%nspinor
    nbnd = header%nband
    nelec = header%nelec
    nspin = nspinor
    IF (nspinor==4) nspin=2
    !
    ncolum = nspin-int(nspinor/nspin)+1
    ALLOCATE( energy(nbnd, ncolum), filling(nbnd, ncolum) )
    !
    ! Make kpoint sets
    ALLOCATE ( link_kpt( nproc, ceiling(nkpt/real(nproc)) ) )
    link_kpt = kpt_array(nkpt, nproc)
    IF (ionode) WRITE(stdout,'(A,I6)') 'Total kptset is ', size(link_kpt,2)
    !
    ! Loop over kpoint sets
    DO ikptset=1,size(link_kpt,2)
      !
      IF (ionode) WRITE(stdout,'(A,I6)') 'kpt set = ', ikptset
      !
      kpt_ind = link_kpt(proc_id+1, ikptset)
      energy = 0d0
      filling = 0d0
      !
      ! skip with kpt index of -1; this proc has no kpt in the set
      IF (kpt_ind/=-1) THEN
        !
        IF (output_type==2 .or. output_type==3) THEN
          IF (allocated(bands)) deallocate(bands)
          IF (allocated(g)) deallocate(g)
        END IF
        !
        ! Read data from the QE output dir
        kdir = kpoint_dir(trim(dir) // '.save', kpt_ind)
        IF (output_type==1 .or. output_type==3) CALL read_eigen_tr(kdir, energy, filling, ncolum)
        IF (output_type==2 .or. output_type==3) CALL read_bands_tr(kdir, bands, nspin)
        IF (output_type==2 .or. output_type==3) CALL read_g_tr(kdir, g)
        !
        IF (output_type==2 .or. output_type==3 ) THEN
          !
          ! Write bands data
          DO ispin=1,nspin
            !
            kpt_file = unk_kpoint_file(trim(dir) // '.unk', kpt_ind, ispin)
            !
            OPEN(iuout, FILE=kpt_file, FORM='unformatted')
            WRITE(iuout) size(bands,1), kpt_ind, nbnd
            DO iband=1,nbnd
              WRITE(iuout) bands(:, ispin, iband)
            END DO
            CLOSE(iuout)
            !
          END DO
          !
          ! Write g vectors
          gvec_file = unk_gvect_file(trim(dir) // '.unk', kpt_ind, 1)
          !
          OPEN(iuout, FILE=gvec_file, FORM='unformatted')
          WRITE(iuout) size(g, 2)
          DO ipw=1,size(g,2)
            write(iuout) g(:,ipw)
          END DO
          CLOSE(iuout)
          !
        END IF
      END IF
      !
      ! Write eigenvalues
      IF (output_type==1 .or. output_type==3) THEN
        !
        ! Gather fillings and energies from all nodes
        IF (ionode) THEN
          IF (allocated(energy_all)) deallocate(energy_all)
          IF (allocated(filling_all)) deallocate(filling_all)
          ALLOCATE ( energy_all( nbnd, ncolum, nproc ), filling_all( nbnd, ncolum, nproc ) )
        END IF
        !
        CALL mp_barrier(mp_comm_all)
        CALL mp_gather(energy, energy_all, ionode_id, mp_comm_all)
        CALL mp_gather(filling, filling_all, ionode_id, mp_comm_all)
        !
        ! writting energy has the unit of Hartree
        IF (ionode) THEN
          DO ispin=1,ncolum
            WRITE(is, fmt='(I1)' ) ispin
            eig_file = trim(dir) // '.unk/' // trim(dir) // '.eigen.' // is
            !
            OPEN(iuout, FILE=eig_file, POSITION='append')
            !
            DO ikpt=1,count(link_kpt(:,ikptset)/=-1)
              DO iband=1,nbnd                 
                WRITE(iuout,'(2I15,2E30.20)') iband, (ikptset-1)*nproc+ikpt, energy_all(iband,ispin,ikpt), filling_all(iband,ispin,ikpt)
              END DO
            END DO
            !
            CLOSE(iuout)
          END DO
        END IF
        !
      END IF
      !
    END DO
    ! END ikptset loop
    !
    IF (ionode) THEN
      WRITE(stdout,'(A)')  '  Writing finished!'
      WRITE(stdout,'(A)') '============================'
    END IF
    !
  END SUBROUTINE trans_qewfc
  !
  ! --------------------------------------------------------------------------
  FUNCTION kpt_array( nkpt, nproc )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp
    !
    IMPLICIT NONE
    !
    INTEGER, ALLOCATABLE :: kpt_array(:,:)
    INTEGER, INTENT(in) :: nkpt, nproc
    !
    INTEGER :: i
    REAL(dp) :: nproc_tmp
    !
    nproc_tmp = real(nproc, kind=dp)
    ALLOCATE ( kpt_array( nproc, ceiling(nkpt/nproc_tmp) ) )
    !
    ! If a processor has no kpt assigned, it gets -1
    kpt_array=-1
    !
    forall(i=1:nkpt) kpt_array(i-(ceiling(i/nproc_tmp)-1)*nproc_tmp, ceiling(i/nproc_tmp))=i
    !
  END FUNCTION kpt_array
  !
  ! --------------------------------------------------------------------------
  FUNCTION kpoint_dir( basedir, ik )
  ! --------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)           :: kpoint_dir
      CHARACTER(LEN=*), INTENT(IN) :: basedir
      INTEGER,          INTENT(IN) :: ik
      !
      CHARACTER(LEN=256) :: kdirname
      CHARACTER(LEN=5)   :: kindex
      CHARACTER(LEN=6)   :: kindex1
      !
      IF (ik<99999) THEN
         WRITE( kindex, FMT = '( I5.5 )' ) ik
         kdirname = TRIM( basedir ) // '/K' // kindex
      ELSEIF (ik<999999) THEN
         WRITE( kindex1, FMT = '( I6.6 )' ) ik
         kdirname = TRIM( basedir ) // '/K' // kindex1
      ELSE
         CALL call_error('kpoint_dir', 'kpoint_dir ik too large, increase format')
      ENDIF
      !
      kpoint_dir = TRIM( kdirname )
      RETURN
      !
  END FUNCTION kpoint_dir
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_eigen_tr( dir, energy, filling, ncolum )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp
    USE io_global, ONLY : iuinput
    USE iotk_module
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: dir
    REAL(dp), INTENT(out) :: energy(:,:), filling(:,:)
    INTEGER, INTENT(in) :: ncolum
    !
    INTEGER :: ispin, nbnd, ierr
    CHARACTER(iotk_attlenx) :: attr
    CHARACTER(LEN=256) :: eigenfiles(3)
    !
    eigenfiles = (/ '//eigenval.xml', '/eigenval1.xml', '/eigenval2.xml'/)
    ierr = 0
    !
    DO ispin=1,ncolum
      !
      CALL iotk_open_read(iuinput, FILE=TRIM(dir)//TRIM(eigenfiles(ispin-1+ncolum)), IERR=ierr)
        CALL iotk_scan_empty(iuinput, "INFO", attr)
        CALL iotk_scan_attr(attr, "nbnd", nbnd)
        CALL iotk_scan_dat(iuinput, "EIGENVALUES", energy(:,ispin) )
        CALL iotk_scan_dat(iuinput, "OCCUPATIONS", filling(:,ispin) )
      CALL iotk_close_read(iuinput)
    END DO
    !
  END SUBROUTINE read_eigen_tr
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_g_tr( dir, g )  
  ! --------------------------------------------------------------------------
    USE io_global, ONLY : iuinput
    USE iotk_module
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: dir
    INTEGER, ALLOCATABLE, INTENT(out)               :: g(:,:)
    !
    INTEGER :: npw
    INTEGER :: ierr
    !
    ierr=0
    !
    CALL iotk_open_read(iuinput, FILE=TRIM(dir)//'/gkvectors.dat', BINARY=.TRUE., IERR=ierr)
    CALL iotk_scan_dat(iuinput, "NUMBER_OF_GK-VECTORS", npw)
    !
    ALLOCATE ( g( 3, npw ) )
    g = 0
    CALL iotk_scan_dat(iuinput, "GRID", g)
    CALL iotk_close_read(iuinput)
    !     
  END SUBROUTINE read_g_tr
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_bands_tr( dir, bands, nspin )
  ! --------------------------------------------------------------------------
    USE io_global, ONLY : iuinput
    USE iotk_module
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: dir
       complex(kind=8), ALLOCATABLE, INTENT(out) :: bands(:,:,:)
    INTEGER, INTENT(in) :: nspin
    !
    INTEGER :: ispin, ibnd, npw, nbnd
    INTEGER :: ierr
    CHARACTER(iotk_attlenx) :: attr
    CHARACTER(LEN=256) :: evcfiles(3)
    !
    evcfiles = (/ '//evc.dat', '/evc1.dat', '/evc2.dat' /)
    ierr = 0
    !
    ! Read dimensions from the data file
    CALL iotk_open_read(iuinput, FILE=TRIM(dir)//TRIM(evcfiles(nspin)), BINARY=.TRUE., IERR=ierr)
      CALL iotk_scan_empty(iuinput, "INFO", attr)
        CALL iotk_scan_attr(attr, "igwx", npw)
        CALL iotk_scan_attr(attr, "nbnd", nbnd)
    CALL iotk_close_read(iuinput)
    !
    ALLOCATE ( bands( npw, nspin, nbnd ) )
    bands = (0.0,0.0)
    !
    ! Read bands from file
    DO ispin=1,nspin
      !
      CALL iotk_open_read(iuinput, FILE=TRIM(dir)//TRIM(evcfiles(ispin-1+nspin)), BINARY=.TRUE., IERR=ierr)
      DO ibnd = 1, nbnd
        CALL iotk_scan_dat(iuinput, "evc" // iotk_index(ibnd), bands(:,ispin,ibnd))
      END DO
      CALL iotk_close_read(iuinput)
      !
    END DO
    !
  END SUBROUTINE read_bands_tr
  !
END MODULE io_wfc

