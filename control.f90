MODULE control
  !
  USE constants, ONLY : dp
  IMPLICIT NONE
  SAVE
  !
  CHARACTER(LEN=256) :: dirname
  CHARACTER(LEN=8) :: data_format
  CHARACTER(LEN=8) :: kpts_type
  !
  REAL(dp) :: cutoff
  ! Photon energy cutoff (Hartree)
  REAL(dp) :: broadwidth
  ! Line broadening (Hartree)
  REAL(dp) :: block_tol
  ! Block degeneracy tolerance (Hartree)
  REAL(dp) :: sv_tol
  REAL(dp) :: sg_tol
  REAL(dp) :: resolution
  !
  LOGICAL :: sg_calc
  ! Do 2nd harmonic gen.
  !
CONTAINS
  !
  SUBROUTINE bcast_control (root)
    !
    USE mp_global, ONLY : mp_bcast, mp_comm_all
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: root
    !
    CALL mp_bcast(dirname, root, mp_comm_all)
    CALL mp_bcast(data_format, root, mp_comm_all)
    CALL mp_bcast(kpts_type, root, mp_comm_all)
    CALL mp_bcast(cutoff, root, mp_comm_all)
    CALL mp_bcast(broadwidth, root, mp_comm_all)
    CALL mp_bcast(block_tol, root, mp_comm_all)
    CALL mp_bcast(sv_tol, root, mp_comm_all)
    CALL mp_bcast(sg_tol, root, mp_comm_all)
    CALL mp_bcast(resolution, root, mp_comm_all)
    CALL mp_bcast(sg_calc, root, mp_comm_all)
    !
  END SUBROUTINE bcast_control
  !
END MODULE control
