MODULE mp_global
  ! Wrappers for MPI usage
  ! --------------------------------------------------------------------------
  USE constants, ONLY : dp
#ifdef __MPI
  USE mpi
#endif
  ! --------------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER, SAVE :: proc_id
  INTEGER, SAVE :: nproc
#ifdef __MPI
  INTEGER, PARAMETER :: mp_comm_all = MPI_COMM_WORLD
#else
  INTEGER, PARAMETER :: mp_comm_all = -1
#endif
  !
  INTERFACE mp_bcast
    MODULE PROCEDURE mp_bcast_int, mp_bcast_ivec, mp_bcast_imat, mp_bcast_it, &
        mp_bcast_rdp, mp_bcast_rvec, mp_bcast_rmat, mp_bcast_rt, &
        mp_bcast_cdp, mp_bcast_cvec, mp_bcast_cmat, mp_bcast_ct, &
        mp_bcast_bool, mp_bcast_char
  END INTERFACE
  !
  INTERFACE mp_reduce_sum
    MODULE PROCEDURE mp_reduce_sum_dp
  END INTERFACE
  !
  INTERFACE mp_gather
    MODULE PROCEDURE mp_gather_vec_dp, mp_gather_mat_dp, mp_gather_dp
  END INTERFACE
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_start
  ! Initialize mpi tasks
  ! --------------------------------------------------------------------------
    USE io_global, ONLY : ionode_id, ionode
    IMPLICIT NONE
#ifdef __MPI
    INTEGER ierr
    !
    CALL MPI_INIT(ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#else
    nproc = 1
    proc_id = 0
#endif
    ionode_id = 0
    ionode = (proc_id == ionode_id)
  END SUBROUTINE mp_start
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_end
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
#ifdef __MPI
    INTEGER ierr
    !
    CALL mp_barrier( MPI_COMM_WORLD )
    CALL MPI_FINALIZE(ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_end
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_stop(errcode)
  ! Abort program with errorcode from MPI wrapper
  ! --------------------------------------------------------------------------
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(in) :: errcode
#ifdef __MPI
    INTEGER :: ierr
#endif
    !
    WRITE(stdout, fmt='( "*** Error in MPI wrapper ***")')
    WRITE(stdout, fmt='( "*** Error code: ",I5)') errcode
#ifdef __MPI
    CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
#endif
    STOP
  END SUBROUTINE mp_stop
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_abort(comm, errcode)
  ! Abort program from elsewhere
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: comm
    INTEGER, INTENT(in) :: errcode
    INTEGER :: ierr
    !
#ifdef __MPI
    CALL MPI_ABORT(comm, errcode, ierr)
#endif
    ERROR STOP
  END SUBROUTINE mp_abort
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_barrier(comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: comm
    INTEGER :: ierr
    !
#ifdef __MPI
    CALL MPI_BARRIER(comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_barrier
  !
  ! --------------------------------------------------------------------------
  !
  ! ... Specific routines for mp_bcast ...
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_bool(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL, INTENT(in) :: buff
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: ierr
    !
    CALL MPI_BCAST(buff, 1, MPI_LOGICAL, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_bool
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_char(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: buff
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: strlen
    INTEGER :: ierr
    !
    strlen = len(buff)
    CALL MPI_BCAST(buff, strlen, MPI_CHARACTER, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_char
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_int(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: buff
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: ierr
    !
    CALL MPI_BCAST(buff, 1, MPI_INTEGER, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_int
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_ivec(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: buff(:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_INTEGER, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_ivec
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_imat(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: buff(:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_INTEGER, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_imat
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_it(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: buff(:,:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_INTEGER, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_it
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_rdp(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: buff
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: ierr
    !
    CALL MPI_BCAST(buff, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_rdp
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_rvec(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: buff(:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_rvec
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_rmat(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: buff(:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_rmat
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_rt(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: buff(:,:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_rt
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_cdp(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(dp), INTENT(in) :: buff
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: ierr
    !
    CALL MPI_BCAST(buff, 1, MPI_DOUBLE_COMPLEX, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_cdp
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_cvec(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(dp), INTENT(in) :: buff(:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_COMPLEX, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_cvec
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_cmat(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(dp), INTENT(in) :: buff(:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_COMPLEX, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_cmat
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_bcast_ct(buff, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(dp), INTENT(in) :: buff(:,:,:)
    INTEGER, INTENT(in) :: root
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    msglen = size(buff)
    CALL MPI_BCAST(buff, msglen, MPI_DOUBLE_COMPLEX, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_bcast_ct
  !
  ! --------------------------------------------------------------------------
  !
  ! ... Specific routines for mp_reduce
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_reduce_sum_dp(buff, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(inout) :: buff(:)
    INTEGER, INTENT(in) :: comm
#ifdef __MPI
    INTEGER :: msglen
    INTEGER :: ierr
    !
    msglen = size(buff)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, buff, msglen, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_reduce_sum_dp
  !
  ! --------------------------------------------------------------------------
  !
  ! ... Specific routines for mp_gather
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_gather_dp(sendbuf, recvbuf, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: sendbuf, recvbuf(:)
    INTEGER, INTENT(in) :: root, comm
#ifdef __MPI
    INTEGER :: msglen, ierr
    !
    msglen = 1
    CALL MPI_GATHER(sendbuf, msglen, MPI_DOUBLE_PRECISION, &
                    recvbuf, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_gather_dp
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_gather_vec_dp(sendbuf, recvbuf, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: sendbuf(:), recvbuf(:,:)
    INTEGER, INTENT(in) :: root, comm
#ifdef __MPI
    INTEGER :: msglen, ierr
    !
    msglen = size(sendbuf)
    CALL MPI_GATHER(sendbuf, msglen, MPI_DOUBLE_PRECISION, &
                    recvbuf, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_gather_vec_dp
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE mp_gather_mat_dp(sendbuf, recvbuf, root, comm)
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: sendbuf(:,:), recvbuf(:,:,:)
    INTEGER, INTENT(in) :: root, comm
#ifdef __MPI
    INTEGER :: msglen, ierr
    !
    msglen = size(sendbuf)
    CALL MPI_GATHER(sendbuf, msglen, MPI_DOUBLE_PRECISION, &
                    recvbuf, msglen, MPI_DOUBLE_PRECISION, root, comm, ierr)
    IF (ierr /= 0) CALL mp_stop( 400 )
#endif
  END SUBROUTINE mp_gather_mat_dp
  !
END MODULE mp_global
