MODULE io_global
  ! Base for file operations
  ! --------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : stdin => INPUT_UNIT, stdout => OUTPUT_UNIT, stderr => ERROR_UNIT
  !
  IMPLICIT NONE
  PRIVATE
  SAVE
  PUBLIC :: stdin, stdout, stderr
  !
  ! ... Set during mp init ...
  PUBLIC :: ionode, ionode_id
  INTEGER :: ionode_id
  LOGICAL :: ionode
  !
  ! --------------------------------------------------------------------------
  ! Units for reading/writitng
  PUBLIC :: iuinput, iuout
  INTEGER :: iuinput = 10
  INTEGER :: iuout = 11
END MODULE io_global
