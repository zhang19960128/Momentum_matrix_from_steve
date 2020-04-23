program pmtxele_program
      use mpi
      use iotk_module
      implicit none

      character(len=256) :: seedname,kdirname
      integer :: bandstart, bandstop

      integer :: iunit, ierr, nkpts, nspin, nbnd, ik,npw,ispin, iband1, iband2, proc_id, nproc, ii, loop_i, loop_j
      logical :: isdatafile

      character(len=5) :: kindex
      character(len=6) :: kindex1

      character(iotk_attlenx)                                                   :: attr
      character(len=256),dimension(3)               :: evcfiles=(/ '/evc.dat', '/evc1.dat', '/evc2.dat' /)

      real(kind=8)                                  :: alat
      real(kind=8), dimension(3,3)                  :: bvecs
      integer, dimension(:,:), allocatable          :: g
      real(kind=8), dimension(:,:), allocatable          :: kpts_cart
      complex(kind=8), dimension(:,:,:), allocatable     :: bands

      real(kind=8), dimension(:,:),  allocatable :: gkvecs

      complex(kind=8),dimension(3) :: pmtxele

      !!Initialize MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,proc_id,ierr)  !!!Absolute rank - myoptid
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)    !!!total number of processes

      namelist /input/ seedname, bandstart, bandstop

      !read arguments
      open(unit=10,file="input_pmat",status="old")
      read(10,nml=input)
      close(10)

      ierr = 0
      iunit = 8

      ! read info
      INQUIRE(FILE=trim(seedname)//'.save/'//'data-file.xml',exist=isdatafile)
      if ( .not. isdatafile) then
        write(*,*) "data-file.xml not found"
      end if
      CALL iotk_open_read( iunit, FILE =  trim(seedname)//'.save/'//'data-file.xml', IERR = ierr )

      if ( ierr > 0 ) write(*,*) "error reading data file"


      CALL iotk_scan_begin( iunit, "BAND_STRUCTURE_INFO" )
      CALL iotk_scan_dat( iunit, "NUMBER_OF_SPIN_COMPONENTS", nspin)
      CALL iotk_scan_dat( iunit, "NUMBER_OF_BANDS", nbnd)
      CALL iotk_scan_end( iunit, "BAND_STRUCTURE_INFO" )

      if (nspin .eq. 4) nspin=2 ! noncollinear spin

      write(*,*) "read nspin nbnd done:" , nspin, nbnd

      !read kpoints
      CALL iotk_scan_begin( iunit, "BRILLOUIN_ZONE" )
      CALL iotk_scan_dat( iunit, "NUMBER_OF_K-POINTS", nkpts)
      allocate(kpts_cart(3,nkpts))
      do ik=1,nkpts
        CALL iotk_scan_empty( iunit, "K-POINT"//TRIM(iotk_index(ik)),   attr )
        CALL iotk_scan_attr(attr, "XYZ",kpts_cart(:,ik)) !2pi/a
      enddo
      CALL iotk_scan_end( iunit, "BRILLOUIN_ZONE" )

      write(*,*) "read kpts done"

      ! read recip latt, a
      CALL iotk_scan_begin( iunit, "CELL" )
      CALL iotk_scan_dat(   iunit, "LATTICE_PARAMETER", alat ) !bohr
      CALL iotk_scan_begin( iunit, "RECIPROCAL_LATTICE_VECTORS" ) !2pi/a
      CALL iotk_scan_dat(   iunit, "b1", bvecs(:,1) )
      CALL iotk_scan_dat(   iunit, "b2", bvecs(:,2) )
      CALL iotk_scan_dat(   iunit, "b3", bvecs(:,3) )
      CALL iotk_scan_end(   iunit, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_scan_end( iunit, "CELL" )

      write(*,*) "read recip done"

      CALL iotk_close_read( iunit )



      open(22,file="pmat.dat")
      !loop over k points
      do ik=1,nkpts
        if (ik<99999) then
          write( kindex, fmt = '( i5.5 )' ) ik
          kdirname = trim( seedname) // '.save/K' // kindex
        elseif (ik<999999) then
          write( kindex1, fmt = '( i6.6 )' ) ik
          kdirname = trim( seedname ) // '.save/K' // kindex1
        else
          write(6,*) 'kpoint_dir ik too large, increase format'
        endif

        !read G vectors
        iunit=100
        ierr=0
        CALL iotk_open_read( iunit, FILE =  TRIM(kdirname)//'/gkvectors.dat', BINARY = .TRUE., IERR = ierr )
        call iotk_scan_dat(   iunit, "NUMBER_OF_GK-VECTORS", npw)
        allocate(g(3,npw))
        g=0
        CALL iotk_scan_dat(   iunit, "GRID", g)
        CALL iotk_close_read( iunit )

        allocate(gkvecs(3,npw))
        gkvecs = spread(kpts_cart(:,ik),2,npw) + g(:,:)
        write(*,*) kpts_cart(:,ik)
        write(*,*) "read gvecs done"
        !debug purpose:!
      !  do loop_i=1,3
      !      do loop_j=1,npw
      !      gkvecs(loop_i,loop_j)=1.0;
      !      END DO
      !  END DO
        !read wavefunctions

        ierr=0
        iunit = 100

        CALL iotk_open_read( iunit, FILE = TRIM(kdirname)//TRIM(evcfiles(nspin)), BINARY = .TRUE., IERR = ierr )
        call iotk_scan_empty(iunit, "INFO", attr)
        call iotk_scan_attr(attr, "igwx",npw)
        call iotk_scan_attr(attr, "nbnd",nbnd)
        CALL iotk_close_read( iunit )

        allocate(bands(npw,nspin,nbnd))
        bands=(0.0,0.0)

        do ispin=1,nspin

          CALL iotk_open_read( iunit, FILE = TRIM(kdirname)//TRIM(evcfiles(ispin-1+nspin)), BINARY = .TRUE., IERR = ierr )
          DO iband1 = 1, nbnd
            CALL iotk_scan_dat( iunit, "evc" // iotk_index( iband1 ), bands(:,ispin,iband1))
          END DO
          CALL iotk_close_read( iunit )

        end do

        write(*,*) "read wfns done"

        !calculate momentum matrix elements
        do iband1=bandstart,bandstop
          do iband2=iband1,bandstop

            !if (nspin.eq.1) then
            !  pmtxele = matmul(gkvecs,conjg(bands(:,1,iband1))*bands(:,1,iband2))
            !else
            pmtxele = matmul(gkvecs,sum(conjg(bands(:,:,iband1))*bands(:,:,iband2),2))
            !endif
            write (22, '(I6, I5, I5, 6(E20.10))' ) ik, iband1, iband2, real(pmtxele), aimag(pmtxele)
          end do
        end do

        write(*,*) "calc p matx done"

        deallocate(g)
        deallocate(gkvecs)
        deallocate(bands)

      end do

      close(22)
end program
