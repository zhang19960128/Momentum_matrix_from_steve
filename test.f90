PROGRAM test_sc
  USE mp_global, ONLY : mp_start, mp_end
  USE control, ONLY: cutoff, block_tol, broadwidth, resolution, sg_calc, dirname
  USE sc_readin, ONLY : read_input
  IMPLICIT NONE
  !
  LOGICAL exst
  ! --------------------------------------------------------------------------
  !
  ! ... Initialize MPI
  CALL mp_start
  !
  ! ... Read control parameters from input
  INQUIRE(FILE='input_sc', EXIST=exst)
  IF (exst) THEN
    CALL read_input
  ELSE
    !!Default parameters
    cutoff = 1.01d0
    ! photon energy cutoff, Hartree
    block_tol = 0.001d0
    ! Block tolerance
    broadwidth = 0.003
    ! spectrum line broadening
    sg_calc = .false.
    ! do 2nd harmonic generation
    resolution = 2.d0
    dirname = "./"
  END IF
  !
  ! Run the sc calculation
  CALL sc_test
  !
  CALL mp_end
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE sc_test
  ! --------------------------------------------------------------------------
  USE constants, ONLY : dp
  USE control, ONLY : cutoff, broadwidth, resolution, sg_calc, dirname
  USE mp_global, ONLY : mp_barrier, mp_comm_all, proc_id
  USE io_global, ONLY : ionode, stdout
  USE es_tools, ONLY : es_header, release_kpt_kl
  USE op_tools, ONLY : kpt_ops, op_blocks_kl, release_matrices_data
  USE sc_tools, ONLY : shift_current_data, shift_kpt, dump_sg_kl, dump_sc_kl
  USE sc_tools, ONLY : energy_grid_ev, release_sc_data, &
     dielectric_spectrum_kl, collect_spectrum_di, write_di_spectrum_kl, &
     shift_current_spectrum2_kl, collect_spectrum_sc, collect_spectrum_sc2, write_sc_spectrum_kl, &
     shift_vector_spectrum_kl, collect_spectrum_sv, write_sv_spectrum_kl, &
     calculate_glass_spectrum, calculate_glass_spectrum2, write_glass_spectrum_kl, &
     sg_spectrum_2omega_kl, collect_spectrum_sg, write_sg_spectrum_kl
  ! --------------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER :: ikptset, ispin
  TYPE(es_header), POINTER :: es_test
  TYPE(kpt_ops) :: matrices
  TYPE(shift_current_data) :: sc
  !
  REAL(dp), ALLOCATABLE :: sc_spec(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: sg_spec(:,:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: pre_sc_spec(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: di_spec(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: sv_spec(:,:,:)
  REAL(dp), ALLOCATABLE :: egrid(:)
  ! --------------------------------------------------------------------------
  !
  es_test=>create_es_test()
  WRITE(stdout,*) proc_id, 'finish creating es_data '
  WRITE(stdout,*) size(es_test%kpt_id_link,1),  size(es_test%kpt_id_link,2), size(es_test%kpt_id_link,3)
  !
  DO ikptset = 1, size(es_test%kpt_id_link,3)
    DO ispin = 1, es_test%nspin
      IF (ionode) WRITE(stdout,*) "kpt_data created for ikptset", ikptset, "ispin=",  ispin
      !
      ! Create test data
      CALL create_wfc_test(es_test, ikptset)
      CALL mp_barrier(mp_comm_all)
      !
      matrices = op_blocks_kl(es_test, 1d0, .true., dirname)
      IF (ionode) WRITE(stdout,*) "matrices created 'op_blocks' for ikptset", ikptset, "ispin=",  ispin
      CALL release_kpt_kl(es_test)
      !
      sc = shift_kpt(matrices, 0.5d0, 0.5d0, es_test%lattice, .true.)
      CALL mp_barrier(mp_comm_all)
      CALL release_matrices_data(matrices)
      !
      ! Prepare spectra
      CALL dielectric_spectrum_kl(sc,di_spec,int(es_test%nspin/2.0)+ispin-1,cutoff,broadwidth,resolution)
      CALL shift_current_spectrum2_kl(sc,egrid,sc_spec,pre_sc_spec,int(es_test%nspin/2.0)+ispin-1,cutoff,broadwidth,resolution)
      CALL shift_vector_spectrum_kl(sc,sv_spec,int(es_test%nspin/2.0)+ispin-1,cutoff,broadwidth,resolution)
      !
      IF (sg_calc==.true.) CALL sg_spectrum_2omega_kl(sc,sg_spec,int(es_test%nspin/2d0)+ispin-1,cutoff,broadwidth,resolution)
      IF (ionode) WRITE(stdout,*) 'written to spectrum for ikptset', ikptset, 'ispin=', ispin
      IF (allocated(sc%block) .and.sg_calc==.false.) CALL dump_sc_kl(sc,int(es_test%nspin/2.0)+ispin)
      IF (allocated(sc%sg_block).and.sg_calc==.true.) CALL dump_sg_kl(sc,int(es_test%nspin/2.0)+ispin)
      !
      CALL release_sc_data(sc)
      IF (ionode) WRITE(stdout,*) 'deallocated sc_date for ikptset', ikptset, 'ispin=', ispin
      IF (ionode) WRITE(stdout,*)
    END DO
    CALL mp_barrier(mp_comm_all)
  END DO
  !
  ! Finalize spectra
  CALL collect_spectrum_di(di_spec)
  CALL collect_spectrum_sv(sv_spec)
  CALL collect_spectrum_sc(sc_spec)
  IF (sg_calc==.true.) CALL collect_spectrum_sg(sg_spec)
  CALL collect_spectrum_sc2(pre_sc_spec)
  !
  ! Write spectra
  IF (ionode) WRITE(stdout,*) "finish calculation and prepare to print out"
  CALL energy_grid_ev(egrid,cutoff,broadwidth,resolution)
  IF (ionode) THEN
    CALL write_sc_spectrum_kl(egrid,sc_spec,di_spec)
    CALL write_di_spectrum_kl(egrid,di_spec)
    CALL write_sv_spectrum_kl(egrid,sv_spec)
    CALL write_glass_spectrum_kl(egrid,sc_spec,di_spec)
    CALL calculate_glass_spectrum2(egrid,pre_sc_spec,di_spec,cutoff,broadwidth,resolution)
    CALL calculate_glass_spectrum(egrid,pre_sc_spec,di_spec,cutoff,broadwidth,resolution)
    IF (sg_calc==.true.) CALL write_sg_spectrum_kl(egrid,sg_spec)
  ENDIF
  !
  END SUBROUTINE sc_test
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE create_wfc_test(es_kl, ikptset)
  ! --------------------------------------------------------------------------
  USE constants, ONLY : dp, dpc, Pi
  USE control, ONLY: cutoff, block_tol, sg_calc
  USE mp_global, ONLY : proc_id
  USE es_tools, ONLY : check_kpt_id_list, kpt_init, es_header
  ! --------------------------------------------------------------------------
  IMPLICIT NONE
  !
  TYPE(es_header), POINTER, INTENT(inout) :: es_kl
  INTEGER, INTENT(in) :: ikptset
  !
  INTEGER :: gvec(3,3)
  REAL(dp), dimension(2) :: filling,energy
  COMPLEX(dpc) :: coeff(3,1,2)
  REAL(dp) :: kcoord(3)
  INTEGER :: kpt_id(2)

  integer :: i, j, ikpt
  ! --------------------------------------------------------------------------
  !
  gvec(:,1) = (/-1, -1, -1/)
  gvec(:,2) = (/0, 0, 0/)
  gvec(:,3) = (/1, 1, 1/)
  energy = (/0d0, 1d0 /)
  filling = (/2d0, 0d0/)
  !
  DO i = -3, 3
    kpt_id = check_kpt_id_list(es_kl%kpt_id_list, es_kl%kpt_id_link(i,proc_id+1,ikptset))
    ikpt = kpt_id(1)
    kcoord = es_kl%kptns(:,1,es_kl%kpt_id_link(i,proc_id+1,ikptset))
    !
    IF (i==0) WRITE(63+proc_id,'(I,I,3F20.10)') proc_id, es_kl%kpt_id_link(0,proc_id+1,ikptset), kcoord
!   coeff(:,1,2)=(/ dcmplx( (1d0/6d0)**0.5,0), dcmplx(0d0,-(1d0/6d0)**0.5), dcmplx(0,(3d0/6d0)**0.5)*exp(dcmplx(0,sum(kcoord)*2*Pi)), dcmplx(0d0,(1d0/6d0)**0.5) /)
!   coeff(:,1,1)=(/ dcmplx(-(3d0/6d0)**0.5,0), dcmplx( (1d0/6d0)**0.5,0d0), dcmplx(0,(1d0/6d0)**0.5)*exp(dcmplx(0,sum(kcoord)*2*Pi)), dcmplx((1d0/6d0)**0.5,0d0) /)
    coeff(:,1,2) = (/ dcmplx( (1d0/3d0)**0.5,0)*exp(dcmplx(0,sum(kcoord)*2*Pi)), &
                      dcmplx( (1d0/3d0)**0.5,0d0), &
                      dcmplx( (1d0/3d0)**0.5,0)*exp(dcmplx(0,sum(kcoord)*2*Pi)) /)
    coeff(:,1,1) = (/ dcmplx( (1d0/6d0)**0.5,0)*exp(dcmplx(0,sum(kcoord)*1*Pi)), &
                      dcmplx(-(4d0/6d0)**0.5,0d0)*exp(dcmplx(0,sum(kcoord)*1*Pi)), &
                      dcmplx((1d0/6d0)**0.5,0)*exp(dcmplx(0,sum(kcoord)*1*Pi)) /)
    CALL kpt_init(kcoord, es_kl%kpt(i), gvec, coeff, energy, filling, block_tol, sg_calc, cutoff)
  END DO
  END SUBROUTINE create_wfc_test
  !
  ! --------------------------------------------------------------------------
  FUNCTION create_es_test
  ! --------------------------------------------------------------------------
  USE constants, ONLY : dp
  USE control, ONLY : cutoff
  USE io_global, ONLY : ionode, stdout
  USE mp_global, ONLY : nproc
  USE es_tools, ONLY : es_header, lattice
  USE es_tools, ONLY : get_kpt_grid_kl, get_kpt_id_list, get_kpt_id_link
  USE flinal_tools, ONLY : ridentity, iidentity
  ! --------------------------------------------------------------------------
  IMPLICIT NONE
  !
  TYPE(es_header), POINTER :: create_es_test
  !
  ! Lattice parameters
  REAL(dp) :: rprimd(3,3)
  INTEGER :: sym(3,3,1)
  REAL(dp) :: lsym(3,1)
  !
  INTEGER :: ikpt, rkpt, skpt, qkpt
  INTEGER :: i, j, k
  ! --------------------------------------------------------------------------
  !
  ALLOCATE(create_es_test)
  !
  ! Create the lattice
  rprimd = 10d0 * ridentity(3)
  sym(:,:,1) = iidentity(3)
  lsym = 0d0
  create_es_test%lattice = lattice(rprimd, sym, lsym)
  !
  create_es_test%freq_cutoff=cutoff
  create_es_test%gen_time_reversal=.false.
  !
  ! Make kpt grid
  create_es_test%nspin = 1
  create_es_test%kptdim = (/8, 8, 8/)
  create_es_test%kptoffset = (/0d0, 0d0, 0d0/)
  create_es_test%nkpt = product(create_es_test%kptdim)
  !
  ALLOCATE( create_es_test%kptns(3,1,create_es_test%nkpt) )
  !
  ikpt = 0
  DO skpt = 1, create_es_test%kptdim(3)
    DO rkpt = 1, create_es_test%kptdim(2)
      DO qkpt = 1, create_es_test%kptdim(1)
        ikpt=ikpt+1
        create_es_test%kptns(:,1,ikpt) = (/ &
            modulo(qkpt+create_es_test%kptdim(1)/2-1,create_es_test%kptdim(1))-create_es_test%kptdim(1)/2+1, &
            modulo(rkpt+create_es_test%kptdim(2)/2-1,create_es_test%kptdim(2))-create_es_test%kptdim(2)/2+1, &
            modulo(skpt+create_es_test%kptdim(3)/2-1,create_es_test%kptdim(3))-create_es_test%kptdim(3)/2+1 /) &
            / (1d0*create_es_test%kptdim)
      END DO
    END DO
  END DO
  !
  CALL get_kpt_grid_kl(create_es_test)
  CALL get_kpt_id_list(create_es_test)
  CALL get_kpt_id_link(create_es_test, nproc)
  !
  IF (ionode) THEN
    WRITE(stdout,*) 'Total Num of kpt set is ', size(create_es_test%kpt_id_link,3), 'with nspin as ', create_es_test%nspin
    WRITE(stdout,*)
    WRITE(60,*) 'nproc= ', nproc, 'nkptset= ', size(create_es_test%kpt_id_link,3)
    WRITE(60,'(7I)') ((create_es_test%kpt_id_link(:,i,j),i=1,nproc),j=1,size(create_es_test%kpt_id_link,3))
    WRITE(61,*) 'nkpt= ', create_es_test%nkpt, 'nsym= ', create_es_test%lattice%nsym, 'Time Rev =', create_es_test%gen_time_reversal
    DO i = 1, create_es_test%nkpt
      WRITE(61,'(50I10)') (create_es_test%kpt_id_list(i,:))
    END DO
    WRITE(62,*) 'nkpt= ', create_es_test%nkpt
    WRITE(62,'(3F20.8)') ((create_es_test%kptns(:,j,i),j=1,create_es_test%lattice%nsym),i=1,create_es_test%nkpt)
  END IF
  END FUNCTION create_es_test
  !
END PROGRAM
