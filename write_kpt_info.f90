
  ! --------------------------------------------------------------------------
  SUBROUTINE write_kpt_info ( header_kl )
  ! --------------------------------------------------------------------------
    USE es_tools, ONLY : es_header
    USE io_global, ONLY : stdout, iuout
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(in) :: header_kl
    CHARACTER(LEN=20) :: fstring
    CHARACTER(LEN=256) :: filename
    INTEGER :: i, j, nproc, nkptset
    !
    WRITE(stdout,*) 'Total Num of kpt set is ', size(header_kl%kpt_id_link,3), 'with nspin as ', header_kl%nspin
    WRITE(stdout,*) 
    CALL flush(stdout)
    !
    filename = "kpt_link.dat"
    nproc = SIZE( header_kl%kpt_id_link, 2 )
    nkptset = SIZE( header_kl%kpt_id_link, 3 )
    OPEN( UNIT = iuout, FILE = trim(filename), FORM = 'FORMATTED', STATUS = 'REPLACE' )
    WRITE(iuout,*) 'nproc= ', nproc, 'nkptset= ', nkptset
    WRITE(iuout,'(7I8)') ((header_kl%kpt_id_link(:,i,j),i=1,nproc),j=1,nkptset)
    CLOSE(iuout, STATUS='keep')
    !
    filename = "kpt_id.dat"
    WRITE(fstring, '("(", I3, "I10)")') header_kl%lattice%nsym
    OPEN( UNIT = iuout, FILE = trim(filename), FORM = 'FORMATTED', STATUS = 'REPLACE' )
    WRITE(iuout,*) 'nkpt= ', header_kl%nkpt, 'nsym= ', header_kl%lattice%nsym, 'Time Rev =', header_kl%gen_time_reversal
    DO i=1,header_kl%nkpt
      WRITE(iuout, fstring) (header_kl%kpt_id_list(i,:))
    END DO
    CLOSE(iuout, STATUS='keep')
    !
    filename = "kpt.dat"
    OPEN( UNIT = iuout, FILE = trim(filename), FORM = 'FORMATTED', STATUS = 'REPLACE' )
    WRITE(iuout,*) 'nkpt= ', header_kl%nkpt
    WRITE(iuout,'(3F20.8)') ((header_kl%kptns(:,j,i),j=1,header_kl%lattice%nsym),i=1,header_kl%nkpt)
    CLOSE(iuout, STATUS='keep')
    !
  END SUBROUTINE write_kpt_info

