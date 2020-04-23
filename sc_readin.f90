MODULE sc_readin
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: read_input
  CHARACTER(LEN=256), SAVE :: input_file = ' '
  !
CONTAINS
  ! --------------------------------------------------------------------------
  SUBROUTINE read_input
  ! --------------------------------------------------------------------------
    !
    USE constants, ONLY : dp, hartree_ev
    USE control, ONLY : dirname, kpts_type, data_format, cutoff, broadwidth, &
                        block_tol, sv_tol, sg_tol, resolution, sg_calc
    USE control, ONLY : bcast_control
    USE io_global, ONLY : stdout, ionode, ionode_id, iuinput
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    NAMELIST /input/ &
        dirname, kpts_type, data_format, cutoff, broadwidth, block_tol, &
        sv_tol, sg_tol, resolution, sg_calc
    ! --------------------------------------------------------------------------
    ! Default parameters
    !
    dirname = ' '
    kpts_type = 'mp' ! monkhorst-pack grid
    data_format = 'auto'
    !
    cutoff = 3.0
    ! photon energy cutoff, Hartree
    broadwidth = 0.003
    ! spectrum line broadening, Hartree
    block_tol = 0.001
    ! Block tolerance, Hartree
    sv_tol = 0.5d0
    sg_tol = 0.5d0
    resolution = 2d0
    sg_calc = .false.
    ! do 2nd harmonic generation 
    !
    ! Read in args
    IF (ionode) THEN
        ierr = open_input_file('input_sc')
        IF (ierr /= 0) CALL call_error("(read_input)", "Error opening input.")
        READ(iuinput, NML=input, IOSTAT=ierr)
        IF (ierr /= 0) CALL call_error("(read_input)", "Error reading input.")
    END IF
    !
    ! Echo output and broadcast to all nodes
    !
    CALL bcast_control(ionode_id)
    !
    IF (ionode) THEN
        WRITE(stdout,*)
        WRITE(stdout,'(A, A   )') "Prefix is                          ", trim(dirname)
        WRITE(stdout,'(A, A   )') "kpts type is                       ", trim(kpts_type)
        WRITE(stdout,'(A,A    )') "Data format is                     ", trim(data_format)
        WRITE(stdout,'(A, F7.4)') "Frequency cutoff (eV)              ", hartree_ev*cutoff
        WRITE(stdout,'(A, F7.4)') "Broadwidth                         ", broadwidth
        WRITE(stdout,'(A, F7.4)') "Degenerate block tolerance         ", block_tol
        WRITE(stdout,'(A, F7.4)') "Shift vector tolerance             ", sv_tol
        WRITE(stdout,'(A, F7.4)') "2nd harmonic gen. tolerance        ", sg_tol
        WRITE(stdout,'(A, F7.4)') "Resolution                         ", resolution
        WRITE(stdout,'(A, L   )') "Calculate SHG response ?           ", sg_calc
        WRITE(stdout,*)
    END IF
  END SUBROUTINE
  !
  !---------------------------------------------------------------------------
  INTEGER FUNCTION open_input_file ( input_file_ )
  ! --------------------------------------------------------------------------
    !
    ! ...  Open file "input_file_" for input read
    ! ...  If "input_file" does not exist, the standard input is dumped to temporary 
    ! ...  file "input_tmp.in"  and this is opened for read
    ! ...  On exit:
    ! ...    Returns -1 if standard input is dumped to file
    ! ...    Returns  0 if input file is successfully opened
    ! ...    Returns  1 if called with wrong arguments
    ! ...    Returns  2 if there was an error opening file
    ! --------------------------------------------------------------------------
    USE io_global, ONLY : stdout, stdin, iuinput
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), intent(in) :: input_file_
    INTEGER :: ierr
    CHARACTER(LEN=512) :: dummy
    LOGICAL :: use_infile
    !
    INQUIRE(FILE=input_file_, EXIST=use_infile)
    !
    IF (use_infile) THEN
      input_file = input_file_
    ELSE
      !
      ! if no file found then copy from standard input
      input_file="input_tmp.in"
      OPEN(UNIT = iuinput, FILE=trim(input_file), FORM='formatted', &
           STATUS='unknown', IOSTAT = ierr )
      IF ( ierr > 0 ) GO TO 30
      !
      dummy=' '
      WRITE(stdout, '(5x,a)') "Waiting for input..."
      DO WHILE ( .true. )
         READ (stdin, fmt='(A512)', IOSTAT=ierr) dummy
         IF (ierr /= 0) EXIT
         WRITE (iuinput,'(A)') trim(dummy)
      END DO
      !
      CLOSE ( UNIT=iuinput, STATUS='keep' )
    END IF
    !
    IF ( input_file .NE. "input_tmp.in") THEN
        WRITE(stdout, '(5x,a)') "Reading input from "//TRIM(input_file)
    ELSE
        WRITE(stdout, '(5x,a)') "Reading input from standard input"
    END IF
    OPEN ( UNIT = iuinput, FILE = TRIM(input_file), FORM = 'FORMATTED', &
           STATUS = 'OLD', IOSTAT = ierr )
    !
    IF ( ierr > 0 ) GO TO 30
    !
    open_input_file = 0
    RETURN
30  open_input_file = 2
    WRITE(stdout, "('Open_input_file: error opening ',A)") TRIM(input_file)
    RETURN
    !
  END FUNCTION open_input_file
  !
  ! --------------------------------------------------------------------------
  INTEGER FUNCTION close_input_file ( )
    !
    ! ...  this subroutine closes the input file opened by open_input_file
    ! ...  removes temporary file if data was read from stdin (text file)
    ! ...  returns -1 if unit is not opened, 0 if no problem, > 0 if problems
    ! --------------------------------------------------------------------------
    !
    USE io_global, ONLY : iuinput
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    INTEGER :: ierr
    !
    INQUIRE ( iuinput, opened = opnd )
    IF (opnd) THEN
      !
      IF ( TRIM(input_file) == "input_tmp.in") THEN
          CLOSE (UNIT=iuinput, STATUS='delete', IOSTAT=ierr )
      ELSE
          CLOSE (UNIT=iuinput, STATUS='keep', IOSTAT=ierr )
      ENDIF
      !
    ELSE
      ierr = -1
    ENDIF 
    !
    close_input_file = ierr
    !
  END FUNCTION close_input_file
  !
END MODULE sc_readin
