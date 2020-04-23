

  SUBROUTINE call_error(arg_1,arg_2)
      USE mp_global, ONLY : mp_abort, mp_comm_all
      USE io_global, ONLY : stdout, ionode
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: arg_1
      CHARACTER(LEN=*), INTENT(in) :: arg_2
integer :: errocode
      !
      IF (ionode) THEN
          WRITE(stdout,*)
          WRITE(stdout,'(A)') repeat('-',48)
          WRITE(stdout,'(A)') "ERROR:" 
          WRITE(stdout,'(A,5X,A)') trim(arg_1), trim(arg_2)
          WRITE(stdout,'(A)') repeat('-',48)
      ENDIF
      !
      CALL mp_abort(mp_comm_all, errocode)
      !
  END SUBROUTINE
  !
  SUBROUTINE call_warning(arg_1,arg_2)
      USE io_global, ONLY : stdout, ionode
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: arg_1
      CHARACTER(LEN=*), INTENT(in) :: arg_2
      !
      IF (ionode) THEN
          WRITE(stdout,*)
          WRITE(stdout,'(A)') repeat('-',48)
          WRITE(stdout,'(A)') "WARNING:" 
          WRITE(stdout,'(A,5X,A)') trim(arg_1), trim(arg_2)
          WRITE(stdout,'(A)') repeat('-',48)
      ENDIF
  END SUBROUTINE

