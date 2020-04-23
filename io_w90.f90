MODULE io_w90
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: read_wakpt, check_unk_files, unk_kpoint_file, unk_gvect_file
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  FUNCTION unk_kpoint_file( basedir, ik ,ispin)
  ! --------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: unk_kpoint_file
    CHARACTER(LEN=*), INTENT(IN) :: basedir
    INTEGER, INTENT(IN) :: ik, ispin
    !
    CHARACTER(LEN=256) :: kfilename
    CHARACTER(LEN=5)   :: kindex
    CHARACTER(LEN=6)   :: kindex1
    CHARACTER(LEN=1)   :: spinindex
    !
    IF (ik<99999) THEN
           WRITE( kindex, FMT = '( I5.5 )' ) ik     
           WRITE( spinindex, FMT = '( I1.1 )' ) ispin     
           kfilename = TRIM( basedir ) // '/UNK' // kindex // '.' // spinindex
    ELSEIF (ik<999999) THEN
           WRITE( kindex1, FMT = '( I6.6 )' ) ik     
           WRITE( spinindex, FMT = '( I1.1 )' ) ispin    
           kfilename = TRIM( basedir ) // '/UNK' // kindex1 // '.' // spinindex
    ELSE
      CALL call_error('unk_kpoint_file', 'kpoint_dir ik too large, increase format')
    ENDIF
    !
    unk_kpoint_file = TRIM( kfilename )
    RETURN
    !
  END FUNCTION unk_kpoint_file
  !
  ! --------------------------------------------------------------------------
  FUNCTION unk_gvect_file( basedir, ik ,ispin)
  ! --------------------------------------------------------------------------
        !
        CHARACTER(LEN=256)           :: unk_gvect_file
        CHARACTER(LEN=*), INTENT(IN) :: basedir
        INTEGER,          INTENT(IN) :: ik, ispin
        !
        CHARACTER(LEN=256) :: gfilename
        CHARACTER(LEN=5)   :: gindex
        CHARACTER(LEN=6)   :: gindex1
        CHARACTER(LEN=1)   :: spinindex
        !
        IF (ik<99999) THEN
           WRITE( gindex, FMT = '( I5.5 )' ) ik     
           WRITE( spinindex, FMT = '( I1.1 )' ) ispin     
           gfilename = TRIM( basedir ) // '/g' // gindex // '.' // spinindex
        ELSEIF (ik<999999) THEN
           WRITE( gindex1, FMT = '( I6.6 )' ) ik     
           WRITE( spinindex, FMT = '( I1.1 )' ) ispin    
           gfilename = TRIM( basedir ) // '/g' // gindex1 // '.' // spinindex
        ELSE
           CALL call_error('unk_gvect_file', 'kpoint_dir ik too large, increase format')
        ENDIF
        !
        unk_gvect_file = TRIM( gfilename )
        !
        RETURN
        !
  END FUNCTION unk_gvect_file
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_wakpt( ikpt, kcoord, fft, energy, filling, nband, ispin, tr)
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp, dpc
    USE control, ONLY : block_tol, cutoff, sg_calc, dirname
    USE es_tools, ONLY : kpt_data, is_tr, kpt_init, kpt_init_tr
    !
    IMPLICIT NONE
    !
    TYPE(kpt_data) :: read_wakpt
    INTEGER, INTENT(in) :: ikpt, ispin, nband
    INTEGER, INTENT(in) :: fft(3)
    REAL(dp), INTENT(in) :: kcoord(3), energy(:), filling(:)
    LOGICAL, INTENT(in) :: tr
    !
    CHARACTER(LEN=256) :: kfile
    COMPLEX(dpc), ALLOCATABLE :: bands(:,:,:)
    INTEGER, ALLOCATABLE :: g(:,:)
    !
    kfile = unk_kpoint_file(trim(dirname) // '.unk', ikpt, ispin)
    CALL read_wannier_bands(kfile, fft, nband, bands, g, ispin, ikpt)
    !
    IF (is_tr(kcoord)) THEN
      CALL kpt_init_tr( kcoord, read_wakpt, g, bands, energy, filling, tr )
    ELSE
      CALL kpt_init( kcoord, read_wakpt, g, bands, energy, filling )
    END IF
    !
    IF (allocated(bands)) deallocate(bands)
    IF (allocated(g)) deallocate(g)
    !
  END FUNCTION read_wakpt
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_wannier_bands( filename, fft, nband, bands, g, ispin, ikpt )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp, dpc
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: ispin, nband, ikpt, fft(3)
    COMPLEX(dpc), ALLOCATABLE, INTENT(out) :: bands(:,:,:)
    INTEGER, ALLOCATABLE, INTENT(out) :: g(:,:)
    !
    COMPLEX(dpc), ALLOCATABLE :: bands_fft(:,:,:,:)
    INTEGER :: iband, ipw, ipwx, ipwy, ipwz, nspinor
    !
    ! Initialize bands and gvecs
    CALL read_bandsfile( filename, fft, ikpt, nband, bands_fft)
    CALL init_wannier_g( fft, g )
    !
    g=-1*g
    DO iband=1,nband
      bands_fft(:,:,:,iband) = fft_to_g( bands_fft(:,:,:,iband) / product(fft) )
    END DO
    !
    !-reorder bands
    !-only for non- spin orbital calc
    !
    nspinor = 1
    ALLOCATE ( bands( product(fft), nspinor, nband ) )
    !
    DO iband=1,nband
      DO ipwz=1,fft(3)
        DO ipwy=1,fft(2)
          DO ipwx=1,fft(1)
            ! Flatten pw index
            ipw = ipwx+(ipwz-1)*fft(3)*fft(2)+(ipwy-1)*fft(2)
            bands(ipw,nspinor,iband) = bands_fft(ipwx,ipwy,ipwz,iband)
          END DO
        END DO
      END DO
    END DO
    !
    IF (allocated(bands_fft)) deallocate(bands_fft)
    !
  END SUBROUTINE read_wannier_bands
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_bandsfile( filename, fft_grid, ikpt, nband, bands )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dpc
    USE io_global, ONLY : iuinput
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(in) :: filename
    INTEGER, INTENT(in) :: ikpt, nband, fft_grid(3)
    COMPLEX(dpc), ALLOCATABLE, INTENT(out) :: bands(:,:,:,:)
    !
    INTEGER :: fft_temp(3)
    INTEGER :: nband_temp, ikpt_temp
    INTEGER :: i, j, k, iband
    !
    ! need to change UNK files to binary.
    OPEN(iuinput, FILE=trim(filename), STATUS='old', FORM='formatted')
    READ(iuinput,*) fft_temp(1), fft_temp(2), fft_temp(3), ikpt_temp, nband_temp
    !
    ! Verify that all indexes and dimensions match
    IF ( nband_temp/=nband .or. any(fft_temp/=fft_grid) .or. ikpt_temp/=ikpt ) THEN
      CALL call_error("read_bandsfile", "nscf and wannier90 not consistent calculation")
    END IF
    !
    ALLOCATE ( bands( fft_grid(1), fft_grid(2), fft_grid(3), nband ) )
    !
    ! the order for reading wfc depends on how it writes (from x->y->z)
    READ(iuinput,*) (((( bands(i,j,k,iband), i=1,fft_grid(1)), j=1,fft_grid(2)), k=1,fft_grid(3)), iband=1,nband )
    CLOSE(iuinput)
    !
  END SUBROUTINE read_bandsfile
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE init_wannier_g( fft_grid, g_vec )
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: fft_grid(3)
    INTEGER, ALLOCATABLE, INTENT(out) :: g_vec(:,:)
    !
    INTEGER :: i, j, k, ipw

    ALLOCATE ( g_vec( 3, product(fft_grid) ) )
    !
    DO k=1,fft_grid(3)
      DO j=1,fft_grid(2)
        DO i=1,fft_grid(1)
          !
          ! Flatten the ijk index (i is fast)
          ipw = i+(j-1)*fft_grid(2)+(k-1)*fft_grid(3)*fft_grid(2)
          g_vec(:,ipw) = grid_vector(fft_grid, (/i,j,k/))
        END DO
      END DO
    END DO
  END SUBROUTINE init_wannier_g
  !
  ! --------------------------------------------------------------------------
  FUNCTION fft_to_g( r_array )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dpc
    !
    IMPLICIT NONE
    !
    COMPLEX(dpc), INTENT(in) :: r_array(:,:,:)
    COMPLEX(dpc), DIMENSION( size(r_array,1), size(r_array,2), size(r_array,3) ) :: fft_to_g
    !
    INTEGER(kind=8) :: plan
    INTEGER :: t_start, t_end
    !
    CALL dfftw_plan_dft_3d(plan,size(r_array,1),size(r_array,2),size(r_array,3),r_array,fft_to_g,1,64)
    CALL dfftw_execute(plan)
    CALL dfftw_destroy_plan(plan)
  END FUNCTION fft_to_g
  !
  ! --------------------------------------------------------------------------
  PURE FUNCTION grid_vector( fft, ng_index )
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: fft(3), ng_index(3)
    INTEGER :: grid_vector(3)
    !
    INTEGER :: ngx,ngy,ngz
    !
    ngz=fft(3)
    ngy=fft(2)
    ngx=fft(1)
    !
    grid_vector=(/ modulo(ng_index(1)+ngx/2-2,ngx)-ngx/2+1, &
                   modulo(ng_index(2)+ngy/2-2,ngy)-ngy/2+1, &
                   modulo(ng_index(3)+ngz/2-2,ngz)-ngz/2+1 /)
  END FUNCTION grid_vector
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE check_unk_files(nkpt,nspin,sporb)
  ! --------------------------------------------------------------------------
    USE control, ONLY : dir => dirname
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkpt, nspin
    LOGICAL, INTENT(in) :: sporb
    !
    INTEGER :: ikpt, ispin, ncolum
    LOGICAL :: kptfile, eigenfile, eigenfile2
    CHARACTER(LEN=256) :: filename
    !
    WRITE(stdout,'(A)') "checking unk files ... "
    !
    ncolum = 1
    IF (sporb) ncolum = 2
    !
    ! Check unk file exist for all spin and kpt
    DO ikpt=1,nkpt
      DO ispin=1,ncolum
        filename = unk_kpoint_file(trim(dir) // '.unk', ikpt, ispin)
        INQUIRE(FILE=filename, EXIST=kptfile)
        IF (.not. kptfile) CALL call_error('check_unk_files', trim(filename) // " is missing !")
      END DO
    END DO
    !
    ! Check eigenfiles exist
    INQUIRE(FILE=trim(dir)//'.unk'//'/'//trim(dir)//'.eigen.1', EXIST=eigenfile)
    IF (nspin==2) THEN
      INQUIRE(FILE=trim(dir)//'.unk'//'/'//trim(dir)//'.eigen.2', EXIST=eigenfile2)
      eigenfile = eigenfile .and. eigenfile2
    END IF
    !
    IF (.not. eigenfile) CALL call_error('check_unk_files', "eigen files is missing !")
    !
    WRITE(stdout,'(A)') "checking unk files complete "
    WRITE(stdout,*)
    !
  END SUBROUTINE check_unk_files
  !
END MODULE io_w90
