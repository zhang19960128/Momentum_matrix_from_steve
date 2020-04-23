
MODULE read_tools
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: header_kl, read_wfc_kl
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  FUNCTION header_kl(output_type)
  ! --------------------------------------------------------------------------
    USE control, ONLY : cutoff, kpts_type, dir => dirname
    USE es_tools, ONLY : es_header, es_header_bcast, energy_filling_bcast, get_kpt_grid_kl, get_kpt_id_list, get_kpt_id_link
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE io_ascii, ONLY : read_es_header_ascii, read_eigens_ascii
    USE io_espresso, ONLY : read_es_header_qe, read_eigens
    USE io_w90, ONLY : check_unk_files
    USE io_wfc, ONLY : trans_qewfc
    USE mp_global, ONLY : nproc, mp_comm_all
    !
    IMPLICIT NONE
    !
    TYPE(es_header), POINTER :: header_kl
    INTEGER, INTENT(inout) :: output_type
    !
    ! ... Read in the header
    IF (ionode) THEN
      IF (output_type == 5) THEN
        header_kl=>read_es_header_ascii()
      ELSE
        header_kl=>read_es_header_qe()
        IF (output_type==2) header_kl%gen_time_reversal = .false.
      END IF
      !
      header_kl%freq_cutoff = cutoff
      header_kl%kpts_type = kpts_type
    END IF
    !
    ! ... Broadcast the header
    CALL es_header_bcast(header_kl, ionode_id, mp_comm_all)
    !
    ! ... If necessary, create the unk files
    IF (output_type == 3) THEN
      CALL trans_qewfc(header_kl, 3)
      output_type = 1
    END IF
    !
    ! ... Read the eigenvalues and fillings
    IF (ionode) THEN
      IF (output_type == 5) THEN
        CALL read_eigens_ascii(header_kl)
      ELSE
        CALL check_unk_files(header_kl%nkpt, header_kl%nspin, header_kl%spinorb)
        CALL read_eigens(dir, header_kl, output_type)
      END IF
    END IF
    !
    ! ... Broadcast the eigenvalues and fillings
    CALL energy_filling_bcast(header_kl, ionode_id, mp_comm_all)
    !
    ! ... Set up the kpt groups
    IF (trim(kpts_type) .ne. 'path') CALL get_kpt_grid_kl(header_kl)
    CALL get_kpt_id_list(header_kl)
    CALL get_kpt_id_link(header_kl, nproc)
    !
    IF (ionode) THEN
      CALL write_kpt_info(header_kl)
    END IF
    !
  END FUNCTION header_kl
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_wfc_kl ( es_kl, ikptset, ispin, output_type )
  ! --------------------------------------------------------------------------
    USE es_tools, ONLY : es_header, check_kpt_id_list, release_kpt_kl
    USE io_global, ONLY : stdout
    USE mp_global, ONLY : proc_id
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(inout) :: es_kl
    INTEGER, INTENT(in) :: ispin, ikptset, output_type
    !
    INTEGER :: kpt_id(2)
    INTEGER :: i, ikpt
    !
    IF (ALLOCATED( es_kl%kpt(0)%block )) CALL release_kpt_kl(es_kl)
    IF (es_kl%kpt_id_link(0,proc_id+1,ikptset)==-1) RETURN
    !
    IF ( trim(es_kl%kpts_type)=='path' ) THEN
      kpt_id = check_kpt_id_list( es_kl%kpt_id_list, es_kl%kpt_id_link(0,proc_id+1,ikptset) )
      ikpt = kpt_id(1)
      !
      es_kl%kpt(0) = read_wfc( kpt_id, es_kl%kptns(:,1,ikpt), es_kl%lattice,es_kl%fftdim, ispin, &
                               es_kl%energy(:,ispin,ikpt), es_kl%filling(:,ispin,ikpt), &
                               es_kl%gen_time_reversal, es_kl%spinorb, output_type )
      !
      WRITE(stdout, 101) ikptset, ispin, proc_id+1, kpt_id, es_kl%kptns(:,1,kpt_id(1)), es_kl%kpt(0)%g%k, es_kl%kpt_id_link(:,proc_id+1,ikptset)
      CALL flush(stdout)
    ELSE
      DO i=-3,3
        !
        kpt_id = check_kpt_id_list( es_kl%kpt_id_list, es_kl%kpt_id_link(i,proc_id+1,ikptset) )
        ikpt = kpt_id(1)
        !
        es_kl%kpt(i) = read_wfc( kpt_id, es_kl%kptns(:,1,ikpt), es_kl%lattice,es_kl%fftdim, ispin, &
                                 es_kl%energy(:,ispin,ikpt), es_kl%filling(:,ispin,ikpt), &
                                 es_kl%gen_time_reversal, es_kl%spinorb, output_type )
        !
        es_kl%kpt(0)%kdiff(i)=es_kl%kpt_id_link(i,proc_id+1,ikptset)
        !
        IF (i==0) THEN
          WRITE(stdout,101) ikptset, ispin, proc_id+1, kpt_id, es_kl%kptns(:,1,kpt_id(1)), es_kl%kpt(0)%g%k, es_kl%kpt_id_link(:,proc_id+1,ikptset)
          CALL flush(stdout)
        END IF
      END DO
    END IF
    !
  101 FORMAT( I6, I2, X, I3, X, "(", 2I5,")", 3X, "(", 3F9.5, ")", X, "(", 3F9.5, ")", 7I6 )
  END SUBROUTINE read_wfc_kl
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_wfc( kpt_id, kcoord, lattice, fft, ispin, energy, filling, time_re, sporb, output_type)
  ! --------------------------------------------------------------------------
    !
    USE constants, ONLY : dp
    USE es_tools, ONLY : kpt_data, lattice_data, kpt_sym, release_kpt
    USE io_ascii, ONLY : read_askpt
    USE io_espresso, ONLY : read_qekpt
    USE io_w90, ONLY : read_wakpt
    !
    IMPLICIT NONE
    !
    TYPE(kpt_data) :: read_wfc
    TYPE(lattice_data), INTENT(in) :: lattice
    INTEGER, INTENT(in) :: kpt_id(2), fft(3)
    INTEGER, INTENT(in) :: ispin, output_type
    REAL(dp), INTENT(in) :: kcoord(3), energy(:), filling(:)
    LOGICAL, INTENT(in) :: time_re, sporb
    !
    TYPE(kpt_data) :: kpt_temp
    INTEGER :: nband, ikpt, isym
    LOGICAL :: tr
    !
    ikpt = kpt_id(1)
    nband = size(energy)
    !
    IF (output_type==1) THEN
      ! Reading QE input
      IF (kpt_id(2)<=lattice%nsym) THEN
        tr = .false.
        isym = kpt_id(2)
      ELSE
        tr = .true.
        isym = kpt_id(2) - lattice%nsym
      END IF
      !
      IF (allocated(kpt_temp%block)) CALL release_kpt(kpt_temp)
      kpt_temp = read_qekpt(ikpt, kcoord, energy, filling, nband, ispin, time_re, sporb)
      read_wfc = kpt_sym(kpt_temp,lattice%gsym(:,:,isym),lattice%lsym(:,isym),lattice%ssym(:,:,isym),tr)
      !
    ELSE IF  (output_type==2) THEN
      ! Reading W90 input
      read_wfc = read_wakpt(ikpt, kcoord, fft, energy, filling, nband, ispin, time_re)
      !
    ELSE IF (output_type==5) THEN
      ! Reading from ascii input
      read_wfc = read_askpt(ikpt, kcoord, energy, filling, nband, time_re)
      !
    END IF
    !
  END FUNCTION read_wfc
  !
END MODULE
