MODULE io_ascii
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: read_es_header_ascii, read_eigens_ascii, read_askpt
  !
CONTAINS
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_es_header_ascii ()
  ! --------------------------------------------------------------------------
    !
    USE control, ONLY : kpts_type
    USE es_tools, ONLY : es_header, lattice
    USE io_global, ONLY : stdout, iuinput
    !
    IMPLICIT NONE
    !
    TYPE(es_header), POINTER :: read_es_header_ascii
    !
    REAL(kind=8), DIMENSION(3) :: a1, a2, a3 
    INTEGER :: nelec, nbnd
    INTEGER, DIMENSION(3) :: kptdim,kptoff
    !
    REAL(kind=8) :: at(3,3)
    INTEGER :: nsym
    INTEGER, ALLOCATABLE :: sym(:,:,:)
    REAL(kind=8), ALLOCATABLE :: latsym(:,:)
    !
    INTEGER :: nkpts, i, j, ik, ierr

    NAMELIST /systemdata/ a1, a2, a3 , nelec, nbnd, kptdim, kptoff
    !
    WRITE(stdout,'(A)') "Reading es_header from ascii file..."
    WRITE(stdout,*)
    !
    OPEN(unit=iuinput, file='systemdata', status='old', IOSTAT=ierr)
    IF (ierr /= 0) CALL call_error('read_es_header_ascii', 'Error opening systemdata')
    READ(iuinput, nml=systemdata, IOSTAT=ierr)
    IF (ierr /= 0) CALL call_error('read_es_header_ascii', 'Error opening systemdata')
    CLOSE(iuinput)
    !
    at(:,1) = a1
    at(:,2) = a2
    at(:,3) = a3
    !
    nsym = 1
    allocate(sym(3,3,nsym),latsym(3,nsym))
    do i=1,3
      latsym(i,1) = 0.0
      do j=1,3
        sym(i,j,1)= 0.0
        if(i==j) sym(i,j,1)= 1.0
      end do
    end do
    !
    ALLOCATE(read_es_header_ascii)
    !
    read_es_header_ascii%gen_time_reversal = .false.
    read_es_header_ascii%spinorb = .false.
    read_es_header_ascii%nspin=1
    read_es_header_ascii%lattice=lattice(at,sym,latsym)
    read_es_header_ascii%nelec=nelec
    read_es_header_ascii%nband=nbnd
    read_es_header_ascii%kptdim=kptdim
    read_es_header_ascii%fftdim = 0
    !
    ! kcoord use conventional fractional k space coordinates
    ! kptns is the crystal coord
    forall(i=1:3) read_es_header_ascii%kptoffset(i)=kptoff(i)/(2d0*kptdim(i))
    nkpts = kptdim(1)*kptdim(2)*kptdim(3)
    !
    ALLOCATE(read_es_header_ascii%kptns(3,1,nkpts))
    read_es_header_ascii%nkpt=nkpts
    !
    OPEN(unit=iuinput,file='kcrys.dat')
    DO ik=1,nkpts
      READ(iuinput,*) read_es_header_ascii%kptns(:,1,ik)
    END DO
    CLOSE(iuinput)
    !
    ! ECHO out
    WRITE(stdout,*) 
    WRITE(stdout,'(A,3I3 )') 'K-grid          = ', read_es_header_ascii%kptdim
    WRITE(stdout,'(A, I8 )') 'Num of kpoints  = ', read_es_header_ascii%nkpt
    WRITE(stdout,'(A, I3 )') 'Num of sym      = ', read_es_header_ascii%lattice%nsym
    WRITE(stdout,'(A, I5 )') 'Num of elecs    = ', read_es_header_ascii%nelec
    WRITE(stdout,'(A,3I5 )') 'FFT-grid        = ', read_es_header_ascii%fftdim
    WRITE(stdout,'(A, L  )') 'Time Re. Sym.   = ', read_es_header_ascii%gen_time_reversal
    WRITE(stdout,*) 
    WRITE(stdout,*)
    !
  END FUNCTION
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE read_eigens_ascii(header_kl)
  ! --------------------------------------------------------------------------
    !
    USE constants, ONLY : dp, hartree_ev
    USE es_tools, ONLY : es_header
    USE io_global, ONLY : iuinput
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(inout) :: header_kl
    !
    REAL(dp) :: occ
    INTEGER :: ikpt, iband, ispin, nelec, nkpt, nspin, nband, occband
    CHARACTER(len=16) :: ikpt_str
    INTEGER :: ierr
    !
    CALL call_warning("read_eigens_ascii", "fillings are determined by number of electrons! Only for insulators!")
    !
    nkpt=header_kl%nkpt
    nelec=header_kl%nelec
    nspin=header_kl%nspin
    nband=header_kl%nband
    IF (mod(int(nelec),2) /= 0) &
        CALL call_error("read_eigens_ascii", "num of elec is not even; ferromagnetism not implented!")
    !
    ALLOCATE( header_kl%energy(nband,nspin,nkpt), header_kl%filling(nband,nspin,nkpt) )
    !
    DO ikpt=1,nkpt
      WRITE(ikpt_str, fmt='(I0)') ikpt
      OPEN(iuinput,file= 'work' // '/k'//trim(ikpt_str)//'/ens',status='old', IOSTAT=ierr)
      IF (ierr /= 0) CALL call_error('read_eigens_ascii', "error opening " // 'work' // '/k'//trim(ikpt_str)//'/ens' )
      DO iband=1,nband
        READ(iuinput,*, IOSTAT=ierr) header_kl%energy(iband,1,ikpt)
        IF (ierr /= 0) CALL call_error('read_eigens_ascii', "error reading eigenfile" )
      END DO
    END DO
    !
    header_kl%energy = header_kl%energy/hartree_ev
    !
    occ=2d0
    occband=int(nelec/2d0)
    header_kl%filling=0.0
    forall(iband=1:occband,ispin=1:nspin,ikpt=1:nkpt) header_kl%filling(iband,ispin,ikpt) = occ
    !
  END SUBROUTINE
  !
  ! --------------------------------------------------------------------------
  FUNCTION read_askpt( ikpt, kcoord, energy, filling, nband, tr )
  ! --------------------------------------------------------------------------
    USE constants, ONLY : dp, dpc
    USE control, ONLY : block_tol, cutoff, sg_calc, dirname
    USE es_tools, ONLY : kpt_data, kpt_init
    USE io_global, ONLY : iuinput
    !
    IMPLICIT NONE
    !
    TYPE(kpt_data) :: read_askpt
    !
    INTEGER, INTENT(in) :: ikpt,  nband
    REAL(dp), INTENT(in) :: kcoord(3), energy(:), filling(:)
    LOGICAL, INTENT(in) :: tr
    !
    INTEGER :: npw, ipw, ib, io
    INTEGER, ALLOCATABLE :: g(:,:)
    REAL(dp), ALLOCATABLE :: btemp(:)
    COMPLEX(dpc), ALLOCATABLE :: bands(:,:,:)
    CHARACTER(len=10) :: ikpt_str
    !
    ! first count the number of planewaves
    OPEN(iuinput, FILE = 'gvecs.dat', STATUS = 'old', FORM = 'formatted' )
    npw = 0
    DO
      READ(iuinput, *, IOSTAT=io)
      IF (io.ne.0) EXIT
      npw = npw+1
    END DO
    CLOSE(iuinput)
    !
    ! read planewaves
    ALLOCATE( g(3, npw) )
    OPEN(iuinput, FILE = 'gvecs.dat', STATUS = 'old', FORM = 'formatted' )
    DO ipw=1,npw
      READ(iuinput,*) g(:,ipw)
    END DO
    CLOSE(iuinput)
    !
    ! read wavefunctions
    ALLOCATE ( bands(npw, 1, nband), btemp(2*nband) )
    WRITE( ikpt_str, FMT = '( I0 )' ) ikpt
    OPEN(iuinput, FILE= 'work' // '/k'//trim(ikpt_str)//'/evecs', STATUS='old', FORM='formatted' )
    DO ipw=1,npw
      READ(iuinput,*) btemp(:)
      DO ib=1,nband
        bands(ipw,1,ib) = btemp(2*ib-1)+(0.0,1.0)*btemp(2*ib)
      END DO
    END DO
    CLOSE(iuinput)
    !
    CALL kpt_init(kcoord, read_askpt, g, bands, energy, filling)
    !
    IF (allocated(g)) deallocate(g)
    IF (allocated(bands)) deallocate(bands)
    IF (allocated(btemp)) deallocate(btemp)
    !
  END FUNCTION read_askpt
  !
END MODULE io_ascii
