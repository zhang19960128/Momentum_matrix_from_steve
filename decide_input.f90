    
  ! --------------------------------------------------------------------------
  SUBROUTINE decide_input(output_format)
  ! ... Decide the electronic structure output to use
  ! ... 1 : QE output used, W90 files present
  ! ... 2 : W90 files used
  ! ... 3 : QE output used, W90 files absent
  ! ... 5 : read from ascii format
  ! ... others are errors
  ! --------------------------------------------------------------------------
    !
    USE control, ONLY : data_format, dir => dirname
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp_global, ONLY : mp_bcast, mp_comm_all
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: output_format
    LOGICAL :: qe_1, qe_2, unk_1, unk_2, qe_3, qe_4, qe, unk
    !
    IF (ionode) THEN
      !
      IF ( TRIM(data_format) == 'ascii' ) THEN
        output_format = 5
        GOTO 10
      END IF
      !
#ifdef __INTEL_COMPILER
      INQUIRE(DIRECTORY=trim(dir)//'.unk',  EXIST=unk)
      INQUIRE(DIRECTORY=trim(dir)//'.save', EXIST=qe )
#elif defined( __GFORTRAN__ )
      INQUIRE(FILE=trim(dir)//'.unk/.',  EXIST=unk)
      INQUIRE(FILE=trim(dir)//'.save/.', EXIST=qe )
#else
#error Unsupported compiler
#endif
      INQUIRE(FILE=trim(dir)//'.unk'//'/'//"g00001.1", EXIST=qe_1)
      INQUIRE(FILE=trim(dir)//'.unk'//'/'//"g000001.1", EXIST=qe_2)
      INQUIRE(FILE=trim(dir)//'.unk'//'/'//"UNK00001.1", EXIST=unk_1)
      INQUIRE(FILE=trim(dir)//'.unk'//'/'//"UNK000001.1", EXIST=unk_2)
      INQUIRE(FILE=trim(dir)//'.save'//'/'//"K00001"//'/'//'gkvectors.dat', EXIST=qe_3)
      INQUIRE(FILE=trim(dir)//'.save'//'/'//"K000001"//'/'//'gkvectors.dat', EXIST=qe_4)
      !
      IF (unk .and. qe) THEN
        IF ((qe_1 .or. qe_2) .and. (unk_1 .or. unk_2)) THEN
          WRITE(stdout,'(A)') "Quantum Espresso output used ."
          output_format = 1
        ELSE IF ((.not. qe_1) .and. (unk_1 .or. unk_2)) THEN
          WRITE(stdout,'(A)') "Wannier90 output used ."
          output_format = 2
        ELSE
          output_format = 0
          CALL call_error('decide_input', 'Cannot decide input type')
        END IF
        !
      ELSE IF (qe .and. (.not. unk)) THEN
        IF (qe_3 .or. qe_4) THEN
          WRITE(stdout,'(A)') "Quantum Espresso output used ."
          output_format = 3
        ELSE
          output_format = -2
          CALL call_error('decide_input', 'Please check *.save folder and have wf_collect=.true.')
        END IF
        !
      ELSE IF (unk .and. (.not. qe)) THEN
        output_format = -1
        CALL call_error('decide_input', 'save folder missing')
        !
      ELSE IF ((.not. qe ).and.( .not. unk)) THEN
        output_format = 0
        CALL call_error('decide_input', 'Cannot decide input type')
        !
      ELSE
        output_format = 0
        CALL call_error('decide_input', 'Cannot decide input type')
        !
      END IF
    END IF
    !
10  CONTINUE
    CALL mp_bcast(output_format, ionode_id, mp_comm_all)
    !
  END SUBROUTINE
