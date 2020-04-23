module sc_tools
    use flinal_tools
    use es_tools
    use op_tools
    implicit none

  type shift_block_data
    integer, dimension(2,2)                                             :: band_range=0
    real(dp)                                                            :: energy=0d0, filling=0d0
    real(dp), dimension(3)                                              :: vec=0d0
    real(dp), dimension(3,3)                                            :: intensity=0d0
    real(dp), dimension(3,3,3)                                          :: current=0d0
  end type

  type sg_block_data
    integer, dimension(:,:),allocatable                                 :: band_range
    real(dp),dimension(:,:),allocatable                                 :: energy, filling
    real(dp), dimension(:,:,:,:),allocatable                            :: sg_intensity
  end type

  type shift_current_data
    integer                                                             :: nblocks
    real(dp), dimension(3)                                              :: k,kcoord
    real(dp)                                                            :: dv
    type(shift_block_data), dimension(:,:),allocatable                  :: block
    type(sg_block_data), dimension(:,:),allocatable                     :: sg_block
  end type





interface shift_kpt
    module procedure shift_kpt_data
end interface

contains




    function sv_phase(p, s_kdiff, p_kdiff, tol)
        complex(dpc), dimension(:,:,:),                         intent(in)      :: p
        complex(dpc), dimension(:,:),                           intent(in)      :: s_kdiff
        complex(dpc), dimension(:,:,:),                         intent(in)      :: p_kdiff
        real(dp),                                               intent(in)      :: tol


        real(dp), dimension(size(p,1),size(p,1))                                :: sv_phase


        complex(dp), dimension(size(p,1),size(p,1))                             :: vec
        integer                                                                 :: i




        vec=0d0
        do i=1,3
            vec=vec+MATMUL(p(:,:,i), MATMUL(s_kdiff(:,:), tranjg(p_kdiff(:,:,i))))
        enddo

        if (abs(dble(trace(mat_log(  matmul(s_kdiff, tranjg(s_kdiff)) )))) < tol*size(vec,1)/2d0) then
            sv_phase=dimag(mat_log(vec))
!            write(7+proc_id,*)"ok"
        else
            sv_phase=0d0
!            write(7+proc_id,*)"too_big"
        endif

    end function

    subroutine sv(matrices,iblock,jblock,vec_right,vec_left,lat,tol)
        real(dp), dimension(:,:,:), allocatable,                intent(out)     :: vec_right,vec_left


        type (kpt_ops),                                         intent(in)      :: matrices
        integer,                                                intent(in)      :: iblock, jblock
        type(lattice_data),                                     intent(in)      :: lat
        real(dp),                                               intent(in)      :: tol


        integer                                                                 :: q

!!maybe temporary. done purely for readability
        type(sc_op)                                                             :: rmat, lmat

        rmat%p=>matrices%block(iblock,jblock)%p
        rmat%p_dk=>matrices%block(iblock,jblock)%p_dk
        rmat%s_dk=>matrices%block(jblock,jblock)%s_dk
        rmat%band_range=matrices%block(iblock,jblock)%band_range

        lmat%p=>matrices%block(jblock,iblock)%p
        lmat%p_dk=>matrices%block(jblock,iblock)%p_dk
        lmat%s_dk=>matrices%block(iblock,iblock)%s_dk
        lmat%band_range=matrices%block(jblock,iblock)%band_range

!!

        allocate(vec_right(rmat%band_range(2,1)-rmat%band_range(1,1)+1, rmat%band_range(2,1)-rmat%band_range(1,1)+1,3), &
                 vec_left (lmat%band_range(2,1)-lmat%band_range(1,1)+1, lmat%band_range(2,1)-lmat%band_range(1,1)+1,3))
        !!Finite difference.  Both forward and backward.  Divide by dk in reciprocal lattice vectors
        do q=1,3
            vec_right(:,:,q)=.5*        ( sv_phase(transform(rmat%p%v,lat%gprimd),rmat%s_dk(q)%s,transform(rmat%p_dk(q)%v,lat%gprimd),tol) &
                                            -sv_phase(transform(rmat%p%v,lat%gprimd),rmat%s_dk(-q)%s,transform(rmat%p_dk(-q)%v,lat%gprimd),tol))/matrices%dk(q)
            vec_left(:,:,q)=.5*transpose(-sv_phase(transform(lmat%p%v,lat%gprimd),lmat%s_dk(q)%s,transform(lmat%p_dk(q)%v,lat%gprimd),tol) &
                                            +sv_phase(transform(lmat%p%v,lat%gprimd),lmat%s_dk(-q)%s,transform(lmat%p_dk(-q)%v,lat%gprimd),tol))/matrices%dk(q)


!       enddo


        enddo

        nullify(rmat%p,rmat%p_dk,rmat%s_dk,lmat%p,lmat%p_dk,lmat%s_dk)

!       write(7+proc_id,*) iblock, jblock
!       call flush(7+proc_id)

    end subroutine


    pure function intensity(p)
        complex(dpc), dimension(:,:,:),                         intent(in)      :: p

        real(dp), dimension(size(p,1),size(p,1),3,3)                            :: intensity


        integer                                                                 :: i,j

        forall(i=1:3,j=1:3) intensity(:,:,i,j)=dble(matmul(p(:,:,i), tranjg(p(:,:,j))))

    end function


    pure FUNCTION pijk(p1,p2,p3)
        complex(dpc),dimension(:,:,:),                          intent(in)      :: p1,p2,p3

        real(dp),dimension(size(p1,1),size(p3,2),3,3,3)                         :: pijk

        integer                                                                 :: s, t, q

        forall(s=1:3,t=1:3,q=1:3) pijk(:,:,s,t,q)=aimag(matmul(matmul(p1(:,:,s),p2(:,:,t)),p3(:,:,q)))

    end FUNCTION


    !###modified in 4/22/12
    FUNCTION shift_kpt_data(matrices,sv_tol,sg_tol,lat,sg)
        type (shift_current_data)                                               :: shift_kpt_data

        type (kpt_ops),                                         intent(in)      :: matrices
        real(dp),                                               intent(in)      :: sv_tol
        real(dp),                                               intent(in)      :: sg_tol
        type (lattice_data),                                    intent(in)      :: lat
        logical,                                                intent(in)      :: sg

        integer                                                                 :: iblock,jblock

        if (.not. allocated(matrices%block)) return

        shift_kpt_data%nblocks=matrices%nblocks
        shift_kpt_data%k=matrices%k
        ! transform back to cartesion
        shift_kpt_data%kcoord=matmul(lat%gprimd,matrices%k)
        shift_kpt_data%dv=matrices%dk(1)*matrices%dk(2)*matrices%dk(3)

        allocate(shift_kpt_data%block(shift_kpt_data%nblocks,shift_kpt_data%nblocks),shift_kpt_data%sg_block(shift_kpt_data%nblocks,shift_kpt_data%nblocks))
        do iblock=1,matrices%nblocks
            do jblock=1,matrices%nblocks
                 if (associated(matrices%block(iblock,jblock)%p_dk))  then
                     if (abs(matrices%block(iblock,jblock)%filling)>0d0) shift_kpt_data%block(iblock,jblock)=sc_block(matrices, iblock,jblock, lat, sv_tol)
                     if (sg) shift_kpt_data%sg_block(iblock,jblock)=sg_block(matrices, iblock,jblock, lat, sg_tol)
                 end if
            enddo
        enddo

    end FUNCTION




    FUNCTION sg_block(matrices,iblock,jblock,lat,sg_tol)
        type(sg_block_data)                                                     :: sg_block

        type (kpt_ops),                                         intent(in)      :: matrices
        integer,                                                intent(in)      :: iblock, jblock
        type(lattice_data),                                     intent(in)      :: lat
        real(dp),                                               intent(in)      :: sg_tol


        integer                                                                 :: kblock, i, j, q

        allocate(sg_block%energy(2,matrices%nblocks),sg_block%filling(2,matrices%nblocks))
        allocate(sg_block%band_range(6,matrices%nblocks))
        allocate(sg_block%sg_intensity(3,3,3,matrices%nblocks))

        sg_block%energy=0d0
        sg_block%filling=0d0
        sg_block%band_range=0
        sg_block%sg_intensity=0d0

        do kblock=1,matrices%nblocks
            sg_block%energy(:,kblock)=(/ matrices%block(iblock,kblock)%energy, matrices%block(kblock,jblock)%energy /)
            sg_block%filling(:,kblock)=(/ matrices%block(iblock,kblock)%filling, matrices%block(kblock,jblock)%filling /)
            sg_block%band_range(1:2,kblock)=matrices%block(iblock,kblock)%band_range(:,1)
            sg_block%band_range(3:4,kblock)=matrices%block(iblock,kblock)%band_range(:,2)
            sg_block%band_range(5:6,kblock)=matrices%block(kblock,jblock)%band_range(:,2)
            sg_block%sg_intensity(:,:,:,kblock)=sg_intense(matrices,iblock,kblock,jblock,lat)
        end do

    end FUNCTION



    FUNCTION sg_intense(matrices, iblock, kblock,jblock,lat)
        type(kpt_ops),                                          intent(in)      :: matrices
        type(lattice_data),                                     intent(in)      :: lat
        integer,                                                intent(in)      :: iblock, jblock , kblock

        integer                                                                 :: i,j,q
        real(dp),dimension(:,:,:,:,:),allocatable                               :: intense
        complex(dpc),dimension(:,:,:),allocatable                               :: v1, v2
        real(dp),dimension(3,3,3)                                               :: sg_intense


        allocate(v1(size(matrices%block(iblock,kblock)%p%v,1),size(matrices%block(iblock,kblock)%p%v,2),3),         &
                 v2(size(matrices%block(kblock,jblock)%p%v,1),size(matrices%block(kblock,jblock)%p%v,2),3),         &
                 intense(size(matrices%block(iblock,kblock)%p%v,1),size(matrices%block(jblock,iblock)%p%v,2),3,3,3) )

        v1=matrices%block(iblock,kblock)%p%v(:,:,:)
        v2=matrices%block(kblock,jblock)%p%v(:,:,:)

        if (iblock==kblock) then
            forall(q=1:3) v1(:,:,q)=cidentity(size(matrices%block(iblock,kblock)%p%v,1))*dcmplx((matrices%k(q)),0)
        else if (kblock==jblock) then
            forall(q=1:3) v2(:,:,q)=cidentity(size(matrices%block(kblock,jblock)%p%v,1))*dcmplx((matrices%k(q)),0)
        end if

        intense=pijk(v1,v2,matrices%block(jblock,iblock)%p%v)
!       if (iblock==2.and. jblock==1.and.kblock==1) write(80+proc_id,'(I,3F15.10,F20.10)') proc_id,matrices%k, intense(:,:,2,1,3)
!       if (iblock==2.and. jblock==1.and.kblock==2) write(79-proc_id,'(I,3F15.10,F20.10)') proc_id,matrices%k, intense(:,:,2,1,3)
!       write(proc_id+40,'(A,3F20.10)') '----k==', matrices%k
!       write(proc_id+40,*) iblock,kblock,jblock
!       write(proc_id+40,'(A,I,9F20.10)') 'v1',size(v1), v1(1,1,1), v1(1,1,2), v1(1,1,3)
!       write(proc_id+40,'(A,I,9F20.10)') 'v2',size(v1), v2(1,1,1), v2(1,1,2), v2(1,1,3)
!       write(proc_id+40,'(A,I,9F20.10)') 'v ',size(matrices%block(jblock,iblock)%p%v), matrices%block(jblock,iblock)%p%v(1,1,1), matrices%block(jblock,iblock)%p%v(1,1,2), matrices%block(jblock,iblock)%p%v(1,1,3)
!       call flush(proc_id+40)
!       do q=1,3
!           write(proc_id+40,'(9E20.10)') ((trace(intense(:,:,i,j,q)),i=1,3),j=1,3)
!       call flush(proc_id+40)
!       end do
!       write(proc_id+40,*)
!       call flush(proc_id+40)
   !    forall(i=1:3,j=1:3,q=1:3) intense(:,:,i,j,q)=2d0
        intense=transform(intense,lat%gprimd)/lat%vol
        forall(i=1:3,j=1:3,q=1:3) sg_intense(i,j,q)=trace(intense(:,:,i,j,q))
!       do q=1,3
!           write(proc_id+40,'(9E20.10)') (((sg_intense(i,j,q)),i=1,3),j=1,3)
!       call flush(proc_id+40)
!       end do
!       write(proc_id+40,*)'----------------------'
!       call flush(proc_id+40)

        deallocate(v1, v2, intense)

    end FUNCTION




    function sc_block(matrices, iblock, jblock, lat,tol)
        type(shift_block_data)                                                  :: sc_block

        type (kpt_ops),                                         intent(in)      :: matrices
        integer,                                                intent(in)      :: iblock, jblock
        real(dp),                                               intent(in)      :: tol
        type(lattice_data),                                     intent(in)      :: lat


        integer                                                                 :: i,j,q
        real(dp), dimension(:,:,:), allocatable                                 :: vec_right,vec_left
        real(dp), dimension(:,:,:,:), allocatable                               :: int_right,int_left

        logical :: fileexist
        character(LEN=3)                                                        :: procstr


        sc_block%energy=matrices%block(iblock,jblock)%energy
        sc_block%filling=matrices%block(iblock,jblock)%filling
        sc_block%band_range=matrices%block(iblock,jblock)%band_range

        call sv(matrices,iblock,jblock,vec_right,vec_left,lat,tol)


        allocate(int_right(sc_block%band_range(2,1)-sc_block%band_range(1,1)+1, &
                           sc_block%band_range(2,1)-sc_block%band_range(1,1)+1,3,3), &
                 int_left(sc_block%band_range(2,2)-sc_block%band_range(1,2)+1, &
                          sc_block%band_range(2,2)-sc_block%band_range(1,2)+1,3,3))



        int_right=intensity(matrices%block(iblock,jblock)%p%v)
        int_left=intensity(matrices%block(jblock,iblock)%p%v)


        !!transform to cartesian
        vec_right=transform(vec_right, lat%rprimd/2d0/PI)
        vec_left=transform(vec_left, lat%rprimd/2d0/PI)

        ! in kspace
        int_right=transform(int_right, lat%gprimd)
        int_left=transform(int_left, lat%gprimd)



        !!!For details see notes

        !!intensity is proportional to the absorption coefficient and is output as such.  Units are applied; division by the photonic flux is implicit
            !!.5 to avoid double counting -  int_right and int_left should have the same trace, but we only  need one
            !! there is no /2d0 in the expression for the dielectric(intensity)
        forall(i=1:3,j=1:3)    sc_block%intensity(i,j)=(h_joule/bohr_meter)**2*sc_block%filling*0.5d0/(lat%vol*bohr_meter**3)*(trace(int_right(:,:,i,j)) + trace(int_left(:,:,i,j)))

        !!add left and right and convert
        forall(q=1:3)   sc_block%vec(q)=(-bohr_meter*(trace(vec_right(:,:,q))+trace(vec_left(:,:,q)))*sign(1d0,sc_block%filling))


        !!multiply intensity times vector and apply necessary units
        forall(i=1:3,j=1:3,q=1:3)  sc_block%current(i,j,q)=sc_block%filling/(lat%vol*bohr_meter**3)*(h_joule/bohr_meter)**2*bohr_meter*(trace(matmul(int_right(:,:,i,j),vec_right(:,:,q)))+trace(matmul(int_left(:,:,i,j),vec_left(:,:,q))))

        !!debug
        !!printing out block matrices for even more detailed information
        !!write(procstr,'(I3.3)')proc_id

        !!if(iblock.eq.44 .and. jblock.eq.51 .and. proc_id.le.10) then
          !!inquire(file="debug"//procstr,exist=fileexist)
          !!if(fileexist) then
            !!open(555,FILE="debug"//procstr,status="old",position="append",action="write")
          !!else
            !!open(555,FILE="debug"//procstr,status="new",action="write")
          !!endif
          !!write(555,*) "kpt:"
          !!write(555,*)  sc_block%filling/(lat%vol*bohr_meter**3)*(h_joule/bohr_meter)**2*bohr_meter*(trace(matmul(int_right(:,:,2,2),vec_right(:,:,2)))+trace(matmul(int_left(:,:,2,2),vec_left(:,:,2))))
          !!write(555,*) "int_right:"
          !!write(555,*) int_right
          !!write(555,*) "int_left:"
          !!write(555,*) int_left
          !!write(555,*) "vec_right:"
          !!write(555,*) vec_right
          !!write(555,*) "vec_left:"
          !!write(555,*) vec_left
          !!close(555)
        !!end if





    end function





    SUBROUTINE release_sc_data(sc)
        type(shift_current_data),                               intent(inout)      :: sc

        integer                                                                    :: ikpt, i, j



        if (allocated(sc%block)) deallocate(sc%block)
        if (allocated(sc%sg_block)) then
            do i = 1, size(sc%sg_block,1)
                do j = 1, size(sc%sg_block,2)
                    if (allocated(sc%sg_block(i,j)%energy)) deallocate(sc%sg_block(i,j)%energy)
                    if (allocated(sc%sg_block(i,j)%filling)) deallocate(sc%sg_block(i,j)%filling)
                    if (allocated(sc%sg_block(i,j)%sg_intensity)) deallocate(sc%sg_block(i,j)%sg_intensity)
                    if (allocated(sc%sg_block(i,j)%band_range)) deallocate(sc%sg_block(i,j)%band_range)
                end do
            end do
            deallocate(sc%sg_block)
        end if

    END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!Output functions

  ! --------------------------------------------------------------------------
  SUBROUTINE open_dump_file( ispin, dumpunit )
  ! --------------------------------------------------------------------------
    USE mp_global, ONLY : proc_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ispin
    INTEGER, INTENT(out) :: dumpunit
    !
    CHARACTER(LEN=3) :: procstr
    CHARACTER(LEN=256) :: dump_file(3)
    !
    dump_file(1) = 'dump'
    dump_file(2) = 'dump_up'
    dump_file(3) = 'dump_dw'
    dumpunit = 20+proc_id
    !
    WRITE(procstr,'(I3.3)') proc_id
    OPEN(dumpunit, FILE=TRIM(dump_file(ispin))//procstr, form='formatted', position='append')
    !
  END SUBROUTINE open_dump_file
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE dump_sc( sc, dumpunit )
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(shift_current_data), INTENT(in) :: sc
    INTEGER, INTENT(in) :: dumpunit
    !
    INTEGER :: iblock, jblock
    !
    WRITE(dumpunit, '("kpt = ",3F10.5)') sc%k
    WRITE(dumpunit, '("kcoord = ",3F10.5)') sc%kcoord
    !
    DO iblock=1,sc%nblocks
      DO jblock=1,sc%nblocks
        IF (all(sc%block(iblock,jblock)%band_range/=0)) THEN
          WRITE(dumpunit,'( 6(1X,I0))') sc%block(iblock,jblock)%band_range, iblock, jblock
          WRITE(dumpunit,'( 2E15.6 )') sc%block(iblock,jblock)%filling, sc%block(iblock,jblock)%energy
          WRITE(dumpunit,'( 3E15.6 )') sc%block(iblock,jblock)%vec
          WRITE(dumpunit,'( 9E15.6 )') sc%block(iblock,jblock)%intensity
          WRITE(dumpunit,'( 9E15.6 )') sc%block(iblock,jblock)%current(:,:,1)
          WRITE(dumpunit,'( 9E15.6 )') sc%block(iblock,jblock)%current(:,:,2)
          WRITE(dumpunit,'( 9E15.6 )') sc%block(iblock,jblock)%current(:,:,3)
        END IF
      END DO
    END DO
    !
  END SUBROUTINE dump_sc
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE dump_sc_kl( sc, ispin )
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(shift_current_data), INTENT(in) :: sc
    INTEGER, INTENT(in) :: ispin
    !
    INTEGER :: dumpunit
    !
    CALL open_dump_file( ispin, dumpunit)
    CALL dump_sc( sc, dumpunit )
    !
  END SUBROUTINE dump_sc_kl
  !
  ! TODO: Much duplicated code from dump_sc
  ! --------------------------------------------------------------------------
  SUBROUTINE dump_sg_kl( sc, ispin )
  ! --------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(shift_current_data), INTENT(in) :: sc
    INTEGER, INTENT(in) :: ispin
    !
    INTEGER :: dumpunit, iblock, jblock, kblock
    !
    CALL open_dump_file( ispin, dumpunit )
    !
    WRITE(dumpunit, '("kpt = ",3F10.5)') sc%k
    WRITE(dumpunit, '("kcoord = ",3F10.5)') sc%kcoord
    !
    DO jblock=1,sc%nblocks
      DO iblock=1,sc%nblocks
        IF (all(sc%block(iblock,jblock)%band_range/=0).and.sc%block(iblock,jblock)%energy<0) THEN
          ! Dump sc
          WRITE(dumpunit,'(6I2   )') sc%block(iblock,jblock)%band_range, iblock, jblock
          WRITE(dumpunit,'(2E15.6)') sc%block(iblock,jblock)%filling, -1d0*sc%block(iblock,jblock)%energy
          WRITE(dumpunit,'(9E15.6)') sc%block(iblock,jblock)%intensity
          !
          ! Dump sg
          DO kblock=1,sc%nblocks
            IF (all(sc%sg_block(iblock,jblock)%band_range/=0)) THEN
              WRITE(dumpunit,'(6I2   )') sc%sg_block(iblock,jblock)%band_range(:,kblock)
              WRITE(dumpunit,'(4E15.6)') sc%sg_block(iblock,jblock)%filling(:,kblock),  sc%sg_block(iblock,jblock)%energy(:,kblock)
              WRITE(dumpunit,'( E15.6)') sc%sg_block(iblock,jblock)%sg_intensity(1,2,3,kblock)
            END IF
          END DO
          !
        END IF
      END DO
    END DO
    !
  END SUBROUTINE dump_sg_kl


    subroutine shift_current_spectrum(shift_current,spectrum,cutoff,broadwidth,resolution)
        real(dp), dimension(:,:,:,:), allocatable                               :: spectrum
        type(shift_current_data), dimension(:)                                  :: shift_current
        real(dp),                                               intent(in)      :: cutoff, broadwidth,resolution


        real(dp), parameter                                                     :: units=.5d0*e_charge*(e_charge/e_mass)**2*(h_hartree/hartree_joule**2)

        integer                                                                 :: ikpt,i,j,q,iblock,jblock

        integer                                                                 :: nbins


        nbins=nint(cutoff/broadwidth*resolution)
        allocate(spectrum(nbins,3,3,3))

        spectrum=0d0
        do ikpt=kstart,kend

            do iblock=1,shift_current(ikpt)%nblocks
                do jblock=1,shift_current(ikpt)%nblocks

                    forall (q=1:3,j=1:3,i=1:3,abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0)   spectrum(:,i,j,q)=spectrum(:,i,j,q)+units/shift_current(ikpt)%block(iblock,jblock)%energy**2*shift_current(ikpt)%dv*delta_broaden(shift_current(ikpt)%block(jblock,iblock)%current(i,j,q),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)


                enddo
            enddo

        enddo

        do q=1,3
            do j=1,3
                do i=1,3
                    call collect_spectrum(spectrum(:,i,j,q))
                enddo
            enddo
        enddo

    end subroutine

!///Modify Glass Coefficient
    subroutine shift_current_spectrum2(shift_current,egrid,spectrum,pre_spectrum,cutoff,broadwidth,resolution)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:), allocatable,                  intent(out) :: spectrum,pre_spectrum
        type(shift_current_data), dimension(:)                                  :: shift_current
        real(dp),                                                   intent(in)  :: cutoff, broadwidth,resolution


        real(dp), parameter                                                     :: units=.5d0*e_charge*(e_charge/e_mass)**2*(h_hartree/hartree_joule**2)
        integer                                                                 :: ikpt,i,j,q,iblock,jblock
        integer                                                                 :: nbins

        nbins=nint(cutoff/broadwidth*resolution)
        allocate(spectrum(nbins,3,3,3),pre_spectrum(nbins,3,3,3))

        spectrum=0d0
        pre_spectrum=0d0
        do ikpt=kstart,kend

            do iblock=1,shift_current(ikpt)%nblocks
                do jblock=1,shift_current(ikpt)%nblocks

                    forall (q=1:3,j=1:3,i=1:3,abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0)   spectrum(:,i,j,q)=spectrum(:,i,j,q)+units/shift_current(ikpt)%block(iblock,jblock)%energy**2*shift_current(ikpt)%dv*delta_broaden(shift_current(ikpt)%block(jblock,iblock)%current(i,j,q),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)

                    forall (q=1:3,j=1:3,i=1:3,abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0)   pre_spectrum(:,i,j,q)=pre_spectrum(:,i,j,q)+units/shift_current(ikpt)%block(iblock,jblock)%energy**2*shift_current(ikpt)%dv*delta_funct(shift_current(ikpt)%block(jblock,iblock)%current(i,j,q),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)

                enddo
            enddo

        enddo

        do q=1,3
            do j=1,3
                do i=1,3
                    call collect_spectrum(spectrum(:,i,j,q))
                    call collect_spectrum(pre_spectrum(:,i,j,q))
                enddo
            enddo
        enddo

    end subroutine

    !###modifled in 4/22/12
!///Modify Glass Coefficient
    SUBROUTINE shift_current_spectrum2_kl(shift_current,egrid,spectrum,pre_spectrum,ispin,cutoff,broadwidth,resolution)
        real(dp), dimension(:),                                   intent(in)    :: egrid
        real(dp), dimension(:,:,:,:,:), allocatable,              intent(inout) :: spectrum,pre_spectrum
        type(shift_current_data)                                                :: shift_current
        real(dp),                                                 intent(in)    :: cutoff, broadwidth,resolution
        integer,                                                  intent(in)    :: ispin


        real(dp), parameter                                                     :: units=.5d0*e_charge*(e_charge/e_mass)**2*(h_hartree/hartree_joule**2)
        integer                                                                 :: ikpt,i,j,q,iblock,jblock
        integer                                                                 :: nbins, st

        nbins=nint(cutoff/broadwidth*resolution)

        if (.not. allocated(spectrum)) then
            st=ceiling(ispin/2.0)+1
            allocate(spectrum(nbins,3,3,3,st),pre_spectrum(nbins,3,3,3,st))

            spectrum=0d0
            pre_spectrum=0d0
        end if

        st=int(ispin/2.0)+1
        if (allocated(shift_current%block)) then
            do iblock=1,shift_current%nblocks
                do jblock=1,shift_current%nblocks

                    forall (q=1:3,j=1:3,i=1:3,abs(shift_current%block(iblock,jblock)%energy) > 0d0)   spectrum(:,i,j,q,st)=spectrum(:,i,j,q,st)+units/shift_current%block(iblock,jblock)%energy**2*shift_current%dv*delta_broaden(shift_current%block(jblock,iblock)%current(i,j,q),abs(shift_current%block(iblock,jblock)%energy),broadwidth,nbins,resolution)

                    forall (q=1:3,j=1:3,i=1:3,abs(shift_current%block(iblock,jblock)%energy) > 0d0)   pre_spectrum(:,i,j,q,st)=pre_spectrum(:,i,j,q,st)+units/shift_current%block(iblock,jblock)%energy**2*shift_current%dv*delta_funct(shift_current%block(jblock,iblock)%current(i,j,q),abs(shift_current%block(iblock,jblock)%energy),broadwidth,nbins,resolution)

                enddo
            enddo
        end if


    end SUBROUTINE


    !
    FUNCTION sg_spectrum_omega(shift_current,iblock,jblock,energy,broadwidth,resolution)
        type(shift_current_data)                                                :: shift_current
        real(dp),                                                 intent(in)    :: energy
        real(dp),                                                 intent(in)    :: broadwidth,resolution
        integer,                                                  intent(in)    :: iblock, jblock

        real(dp), dimension(3,3,3,2)                                            :: sg_spectrum_omega

        integer                                                                 :: i, j, q, kblock

!       if (allocated(shift_current%sg_block)) then
            do kblock=1,shift_current%nblocks
                forall (q=1:3,j=1:3,i=1:3)
                    sg_spectrum_omega(i,j,q,1)=sg_spectrum_omega(i,j,q,1)+shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock)                            &
                                              *( ndelta_broaden_value(-1d0*(shift_current%sg_block(iblock,jblock)%energy(1,kblock)),energy,broadwidth,resolution)       &
                                                  *   (shift_current%sg_block(iblock,jblock)%filling(1,kblock))                                                       &
                                               + ndelta_broaden_value(-1d0*(shift_current%sg_block(iblock,jblock)%energy(2,kblock)),energy,broadwidth,resolution)       &
                                                  *   (shift_current%sg_block(iblock,jblock)%filling(2,kblock)) )

                    sg_spectrum_omega(i,j,q,2)=sg_spectrum_omega(i,j,q,2)+shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock)                            &
                                              *(    npp_broaden_value(-1d0*(shift_current%sg_block(iblock,jblock)%energy(1,kblock)),energy,broadwidth,resolution)       &
                                                  *   (shift_current%sg_block(iblock,jblock)%filling(1,kblock))                                                       &
                                               +    npp_broaden_value(-1d0*(shift_current%sg_block(iblock,jblock)%energy(2,kblock)),energy,broadwidth,resolution)       &
                                                  *   (shift_current%sg_block(iblock,jblock)%filling(2,kblock)) )
                end forall
            end do
!       end if

    end FUNCTION



    !
    SUBROUTINE sg_spectrum_2omega(shift_current,spectrum_2omega,ispin,cutoff,broadwidth,resolution)
        type(shift_current_data)                                                :: shift_current
        real(dp),                                                 intent(in)    :: cutoff, broadwidth,resolution
        integer,                                                  intent(in)    :: ispin
        real(dp), dimension(:,:,:,:,:,:), allocatable,            intent(inout) :: spectrum_2omega


        real(dp), dimension(3,3,3,2)                                            :: spectrum_omega

    !   real(dp), parameter                                                     :: units=(e_charge/e_mass)**3*(h_hartree**3)/((2d0*pi)**6)*(h_hartree**3)/(bohr_meter**6)*hartree_joule
        real(dp), parameter                                                     :: units=(e_charge/e_mass)**3/(hartree_joule)**5*(h_joule)**6/((2d0*pi)**3)/(bohr_meter)**6
        real(dp)                                                                :: energy_ij, filling_ij
        integer                                                                 :: ikpt,i,j,q,iblock,jblock, nblock
        integer                                                                 :: nbins, st


        nbins=nint(cutoff/broadwidth*resolution)
        nblock=shift_current%nblocks

        if (.not. allocated(spectrum_2omega)) then
            st=ceiling(ispin/2.0)+1
            allocate(spectrum_2omega(nbins,3,3,3,2,st))

            spectrum_2omega=0d0
        end if

        st=int(ispin/2.0)+1
        if (allocated(shift_current%sg_block)) then
            do jblock=1,nblock
                do iblock=1,nblock
                    if (allocated(shift_current%sg_block(iblock,jblock)%energy)) then
                         energy_ij=-shift_current%sg_block(iblock,jblock)%energy(1,1)-shift_current%sg_block(iblock,jblock)%energy(2,1)
                         filling_ij=-shift_current%sg_block(iblock,jblock)%filling(1,2)-shift_current%sg_block(iblock,jblock)%filling(2,2)
                         if ( energy_ij>0d0 ) then
                              spectrum_omega=sg_spectrum_omega(shift_current,iblock,jblock,energy_ij/2d0,broadwidth,resolution)
                              forall(q=1:3,j=1:3,i=1:3)
                                  spectrum_2omega(:,i,j,q,1,st)=spectrum_2omega(:,i,j,q,1,st)+delta_broaden(spectrum_omega(i,j,q,1),energy_ij,broadwidth,nbins,resolution)*shift_current%dv/(energy_ij)*units/energy_ij**2
                                  spectrum_2omega(:,i,j,q,2,st)=spectrum_2omega(:,i,j,q,2,st)+delta_broaden(spectrum_omega(i,j,q,2),energy_ij,broadwidth,nbins,resolution)*shift_current%dv/(energy_ij)*units/energy_ij**2
                              end forall
                         end if
                    end if
                end do
            end do
        end if

    end SUBROUTINE

    SUBROUTINE sg_spectrum_2omega_kl(shift_current,spectrum_2omega,ispin,cutoff,broadwidth,resolution)
        type(shift_current_data)                                                :: shift_current
        real(dp),                                                 intent(in)    :: cutoff, broadwidth,resolution
        integer,                                                  intent(in)    :: ispin
        real(dp), dimension(:,:,:,:,:,:), allocatable,            intent(inout) :: spectrum_2omega

        real(dp), dimension(3,3,3,2)                                            :: spectrum_omega
        real(dp), parameter                                                     :: units=-(e_charge/e_mass)**3/(hartree_joule)**5*(h_joule)**6/((2d0*pi)**3)/(bohr_meter)**6
        real(dp)                                                                :: energy_ij, filling_ij, omega_delta_kj, omega_delta_ik, omega_inver_kj, omega_inver_ik
        integer                                                                 :: ikpt,i,j,q,iblock,jblock, nblock, kblock
        integer                                                                 :: nbins, st

        real(dp),dimension(:),allocatable                                       :: temp

        nbins=nint(cutoff/broadwidth*resolution)
        nblock=shift_current%nblocks

        allocate(temp(nbins))
        if (.not. allocated(spectrum_2omega)) then
            st=ceiling(ispin/2.0)+1
            allocate(spectrum_2omega(nbins,3,3,3,2,st))

            spectrum_2omega=0d0
        end if

        st=int(ispin/2.0)+1
        if (allocated(shift_current%sg_block)) then
            do jblock=1,nblock
                do iblock=1,nblock
                    if (allocated(shift_current%sg_block(iblock,jblock)%energy)) then
                         energy_ij=-shift_current%sg_block(iblock,jblock)%energy(1,1)-shift_current%sg_block(iblock,jblock)%energy(2,1)
                         filling_ij=-shift_current%sg_block(iblock,jblock)%filling(1,2)-shift_current%sg_block(iblock,jblock)%filling(2,2)
                         if ( energy_ij>0d0 ) then
                              do kblock=1,nblock
                                  omega_delta_kj=ndelta_broaden_value(-shift_current%sg_block(iblock,jblock)%energy(2,kblock),energy_ij/2d0,broadwidth,resolution)
                                  omega_delta_ik=ndelta_broaden_value(-shift_current%sg_block(iblock,jblock)%energy(1,kblock),energy_ij/2d0,broadwidth,resolution)
                                  omega_inver_kj=npp_broaden_value(-shift_current%sg_block(iblock,jblock)%energy(2,kblock),energy_ij/2d0,broadwidth,resolution)
                                  omega_inver_ik=npp_broaden_value(-shift_current%sg_block(iblock,jblock)%energy(1,kblock),energy_ij/2d0,broadwidth,resolution)
                                ! write(70+proc_id,'(3F10.5,3I5,4E20.10)') shift_current%k, iblock,jblock,kblock, omega_delta_kj,omega_delta_ik,omega_inver_kj,omega_inver_ik
                                  temp=pp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_inver_kj
                                ! write(6,*) '!!!!', temp(1001)
                                  forall(q=1:3,j=1:3,i=1:3)
                                      spectrum_2omega(:,i,j,q,1,st)=spectrum_2omega(:,i,j,q,1,st)+                                                                  &
                 !                         ( (shift_current%sg_block(iblock,jblock)%filling(2,kblock))*(-delta_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_delta_kj*pi**2+pp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_inver_kj)    &
                 !                         + (shift_current%sg_block(iblock,jblock)%filling(1,kblock))*(-delta_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_delta_ik*pi**2+pp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_inver_ik) )  &
                 !                         ( (shift_current%sg_block(iblock,jblock)%filling(2,kblock))*(-delta_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_delta_kj*pi**2)    &
                 !                         + (shift_current%sg_block(iblock,jblock)%filling(1,kblock))*(-delta_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,broadwidth,nbins,resolution)*omega_delta_ik*pi**2) )  &
                                           ( (shift_current%sg_block(iblock,jblock)%filling(2,kblock))*(-sg_dd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(2,kblock),broadwidth,nbins,resolution)*pi**2+sg_pp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(2,kblock),broadwidth,nbins,resolution))    &
                                           + (shift_current%sg_block(iblock,jblock)%filling(1,kblock))*(-sg_dd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(1,kblock),broadwidth,nbins,resolution)*pi**2+sg_pp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(1,kblock),broadwidth,nbins,resolution)) )  &
                 !                         ( (shift_current%sg_block(iblock,jblock)%filling(2,kblock))*(-sg_dd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(2,kblock),broadwidth,nbins,resolution)*pi**2)    &
                 !                         + (shift_current%sg_block(iblock,jblock)%filling(1,kblock))*(-sg_dd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(1,kblock),broadwidth,nbins,resolution)*pi**2) )  &
                                          *shift_current%dv/(energy_ij)*units/(energy_ij/2d0)**2

                                      spectrum_2omega(:,i,j,q,2,st)=spectrum_2omega(:,i,j,q,2,st)+                                                                  &
                                           ( (shift_current%sg_block(iblock,jblock)%filling(2,kblock))*(-sg_pd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(2,kblock),broadwidth,nbins,resolution)*pi   -sg_dp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(2,kblock),broadwidth,nbins,resolution)*pi)    &
                                           + (shift_current%sg_block(iblock,jblock)%filling(1,kblock))*(-sg_pd_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(1,kblock),broadwidth,nbins,resolution)*pi   -sg_dp_broaden(shift_current%sg_block(iblock,jblock)%sg_intensity(i,j,q,kblock),energy_ij,-shift_current%sg_block(iblock,jblock)%energy(1,kblock),broadwidth,nbins,resolution)*pi) )  &
                                          *shift_current%dv/(energy_ij)*units/(energy_ij/2d0)**2
                                  end forall
                              end do
                         end if
                    end if
                end do
            end do
        end if

    end SUBROUTINE


    pure function sg_dd_broaden(y,freq,freq_1,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq, freq_1, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: sg_dd_broaden
        integer                                                                 :: nwidth
        integer                                                                 :: i, eindex

        nwidth=4*resolution
        sg_dd_broaden=0d0
        eindex=NINT(freq/broadwidth*resolution)
        do i=max(eindex-nwidth,1),min(eindex+nwidth,nbins)
            sg_dd_broaden(i)=y*(.5d0*broadwidth)/(pi*(((eindex-i)/resolution)**2+.25d0)*broadwidth**2)*ndelta_broaden_value(freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution)
        end do
    end function

    pure function sg_pp_broaden(y,freq,freq_1,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq,freq_1, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: sg_pp_broaden
        integer                                                                 :: i, eindex, nwidth

        nwidth=4*resolution
        sg_pp_broaden=0d0
        eindex=NINT(freq/broadwidth*resolution)
        do i=1,nbins
            sg_pp_broaden(i)=y*((((i-eindex)/resolution)*broadwidth)/( (((eindex-i)/resolution)**2+.25d0)*broadwidth**2 ) - (((i+eindex)/resolution)*broadwidth)/( (((eindex+i)/resolution)**2+.25d0)*broadwidth**2 )) * (npp_broaden_value(freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution)-npp_broaden_value(-freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution))
        end do
    end function

    pure function sg_dp_broaden(y,freq,freq_1,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq, freq_1, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: sg_dp_broaden
        integer                                                                 :: nwidth
        integer                                                                 :: i, eindex

        nwidth=4*resolution
        sg_dp_broaden=0d0
        eindex=NINT(freq/broadwidth*resolution)
        do i=1,nbins
            sg_dp_broaden(i)=y*((((i-eindex)/resolution)*broadwidth)/( (((eindex-i)/resolution)**2+.25d0)*broadwidth**2 ) - (((i+eindex)/resolution)*broadwidth)/( (((eindex+i)/resolution)**2+.25d0)*broadwidth**2 )) * ndelta_broaden_value(freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution)
        end do
    end function

    pure function sg_pd_broaden(y,freq,freq_1,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq,freq_1, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: sg_pd_broaden
        integer                                                                 :: i, eindex, nwidth

        nwidth=4*resolution
        sg_pd_broaden=0d0
        eindex=NINT(freq/broadwidth*resolution)
        do i=max(eindex-nwidth,1),min(eindex+nwidth,nbins)
            sg_pd_broaden(i)=y*(.5d0*broadwidth)/(pi*(((eindex-i)/resolution)**2+.25d0)*broadwidth**2)*(npp_broaden_value(freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution)-npp_broaden_value(-freq_1,i/2d0/resolution*broadwidth,broadwidth,resolution))

        end do
    end function





!///Modify Glass Coefficient
    pure function delta_funct(y,freq,broadwidth,nbins,resolution)
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: y, freq, broadwidth, resolution

        real(dp),dimension(nbins)                                               :: delta_funct
        real(dp)                                                                :: energy_space
        integer                                                                 :: i, eindex

        delta_funct=0
        eindex=NINT(freq/broadwidth*resolution)
        delta_funct(eindex)=y

    end function

!///Modify Glass Coefficient
    pure function delta_broaden_index(y,eindex,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, broadwidth
        integer,                                                intent(in)      :: nbins, eindex
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: delta_broaden_index


        integer                                                                 :: nwidth
        integer                                                                 :: i

        nwidth=4*resolution
        delta_broaden_index=0d0

        forall (i=max(eindex-nwidth,1):min(eindex+nwidth,nbins)) delta_broaden_index(i)=y*(.5d0*broadwidth)/(pi*(((eindex-i)/resolution)**2+.25d0)*broadwidth**2)

    end function


    subroutine shift_vector_spectrum(shift_current,spectrum,cutoff,broadwidth,resolution)
        real(dp), dimension(:,:), allocatable                                   :: spectrum
        type(shift_current_data), dimension(:)                                  :: shift_current
        real(dp),                                               intent(in)      :: cutoff, broadwidth, resolution


        real(dp), parameter                                                     :: units=1d10

        integer                                                                 :: ikpt,i,j,q,iblock,jblock
        real(dp)                                                                :: temp

        integer                                                                 :: nbins
        integer, dimension(2)                                                   :: results


        nbins=nint(cutoff/broadwidth*resolution)
        allocate(spectrum(nbins,3))

        spectrum=0d0
        do ikpt=kstart,kend

            do iblock=1,shift_current(ikpt)%nblocks
                do jblock=1,shift_current(ikpt)%nblocks

!                   forall (q=1:3,abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0)
                       if ( abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0 ) then
                    do q = 1, 3
                          spectrum(:,q)=spectrum(:,q)+units*shift_current(ikpt)%dv*delta_broaden(shift_current(ikpt)%block(jblock,iblock)%vec(q),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)
                    end do
!                   temp=maxval(spectrum)
!                   results=maxloc(spectrum)

!                   write(50+proc_id,*) proc_id, ikpt, iblock, jblock, temp, results, shift_current(ikpt)%block(jblock,iblock)%vec(results(2)), units*shift_current(ikpt)%dv, "=" ,delta_broaden(shift_current(ikpt)%block(jblock,iblock)%vec(results(2)),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)
                       end if


                enddo
            enddo

        enddo

        do q=1,3
            call collect_spectrum(spectrum(:,q))
        enddo

    end subroutine

    !###modifled in 4/22/12
    SUBROUTINE shift_vector_spectrum_kl(shift_current,spectrum,ispin,cutoff,broadwidth,resolution)
        real(dp), dimension(:,:,:), allocatable,                intent(inout)   :: spectrum
        type(shift_current_data)                                                :: shift_current
        real(dp),                                               intent(in)      :: cutoff, broadwidth, resolution
        integer,                                                intent(in)      :: ispin


        real(dp), parameter                                                     :: units=1d10

        integer                                                                 :: ikpt,i,j,q,iblock,jblock
        real(dp)                                                                :: temp

        integer                                                                 :: nbins, st
        integer, dimension(2)                                                   :: results


        nbins=nint(cutoff/broadwidth*resolution)
        if (.not. allocated(spectrum)) then
            st=ceiling(ispin/2.0)+1
            allocate(spectrum(nbins,3,st))

            spectrum=0d0
        end if

        st=int(ispin/2.0)+1
        if (allocated(shift_current%block)) then
            do iblock=1,shift_current%nblocks
                do jblock=1,shift_current%nblocks

!                   forall (q=1:3,abs(shift_current(ikpt)%block(iblock,jblock)%energy) > 0d0)
                    if ( abs(shift_current%block(iblock,jblock)%energy) > 0d0 ) then
                       do q = 1, 3
                          spectrum(:,q,st)=spectrum(:,q,st)+units*shift_current%dv*delta_broaden(shift_current%block(jblock,iblock)%vec(q),abs(shift_current%block(iblock,jblock)%energy),broadwidth,nbins,resolution)
                       end do
                    end if


                enddo
            enddo
        end if


    end SUBROUTINE


    subroutine dielectric_spectrum(shift_current,spectrum,cutoff,broadwidth,resolution)
        real(dp), dimension(:,:,:,:), allocatable                               :: spectrum
        type(shift_current_data), dimension(:)                                  :: shift_current
        real(dp),                                               intent(in)      :: cutoff, broadwidth, resolution


        real(dp), parameter                                                     :: units=1/eps_0*((h_hartree/(2d0*pi))**2/hartree_joule)*(e_charge/e_mass)**2

        integer                                                                 :: ikpt,i,j,q,iblock,jblock
        integer                                                                 :: nbins



        nbins=nint(cutoff/broadwidth*resolution)
        allocate(spectrum(2,nbins,3,3))

        spectrum(1,:,:,:)=0d0
        spectrum(2,:,:,:)=0d0
        do ikpt=kstart,kend

            do iblock=1,shift_current(ikpt)%nblocks
                do jblock=1,shift_current(ikpt)%nblocks

                    forall (i=1:3,j=1:3,shift_current(ikpt)%block(jblock,iblock)%filling > 0d0)   spectrum(1,:,i,j)=spectrum(1,:,i,j)+1d0/(shift_current(ikpt)%block(iblock,jblock)%energy**2)*units*shift_current(ikpt)%dv*(pp_broaden(shift_current(ikpt)%block(jblock,iblock)%intensity(i,j),shift_current(ikpt)%block(iblock,jblock)%energy,broadwidth,nbins,resolution)-pp_broaden(shift_current(ikpt)%block(jblock,iblock)%intensity(i,j),-shift_current(ikpt)%block(iblock,jblock)%energy,broadwidth,nbins,resolution))
                    forall (i=1:3,j=1:3,shift_current(ikpt)%block(jblock,iblock)%filling > 0d0)   spectrum(2,:,i,j)=spectrum(2,:,i,j)+pi/(shift_current(ikpt)%block(iblock,jblock)%energy**2)*units*shift_current(ikpt)%dv*delta_broaden(shift_current(ikpt)%block(jblock,iblock)%intensity(i,j),abs(shift_current(ikpt)%block(iblock,jblock)%energy),broadwidth,nbins,resolution)

                enddo
            enddo

        enddo
        do j=1,3
            do i=1,3
                call collect_spectrum(spectrum(1,:,i,j))
                call collect_spectrum(spectrum(2,:,i,j))
            enddo
        enddo

        !!Real part of dielectric needs +1
        spectrum(1,:,i,j)=spectrum(1,:,i,j)+1d0

    end subroutine

    !###modified in 4/22/12
  ! --------------------------------------------------------------------------
  SUBROUTINE dielectric_spectrum_kl(shift_current,spectrum,ispin,cutoff,broadwidth,resolution)
  ! --------------------------------------------------------------------------
    USE constants, ONLY : pi, eps_0, h_hartree, hartree_joule, e_charge, e_mass
    !
        real(dp), dimension(:,:,:,:,:), allocatable,            intent(inout)   :: spectrum
        type(shift_current_data)                                                :: shift_current
        real(dp),                                               intent(in)      :: cutoff, broadwidth, resolution
        integer,                                                intent(in)      :: ispin

    REAL(dp), PARAMETER :: units=1/eps_0*((h_hartree/(2d0*pi))**2/hartree_joule)*(e_charge/e_mass)**2
    INTEGER :: i,j,iblock,jblock
    INTEGER :: nbins, st
    !
        nbins=nint(cutoff/broadwidth*resolution)
        if (.not. allocated(spectrum)) then
            allocate(spectrum(2,nbins,3,3,(ceiling(ispin/2.0)+1)))

            spectrum(1,:,:,:,:)=0d0
            spectrum(2,:,:,:,:)=0d0
        end if
    !
        st=int(ispin/2.0)+1
    IF ( allocated(shift_current%block)) THEN
      DO iblock = 1, shift_current%nblocks
        DO jblock = 1, shift_current%nblocks
          !
          FORALL (i=1:3, j=1:3, shift_current%block(jblock,iblock)%filling > 0d0) &
            spectrum(1,:,i,j,st) = spectrum(1,:,i,j,st) + units*shift_current%dv &
                * 1d0/(shift_current%block(iblock,jblock)%energy**2) &
                * (pp_broaden(shift_current%block(jblock,iblock)%intensity(i,j), &
                              shift_current%block(iblock,jblock)%energy, &
                              broadwidth, nbins, resolution) &
                 - pp_broaden(shift_current%block(jblock,iblock)%intensity(i,j), &
                              -shift_current%block(iblock,jblock)%energy, &
                              broadwidth, nbins, resolution))
          !
          FORALL (i=1:3, j=1:3, shift_current%block(jblock,iblock)%filling > 0d0) &
            spectrum(2,:,i,j,st) = spectrum(2,:,i,j,st) + units*shift_current%dv &
                * pi/(shift_current%block(iblock,jblock)%energy**2) &
                * delta_broaden(shift_current%block(jblock,iblock)%intensity(i,j), &
                                abs(shift_current%block(iblock,jblock)%energy), &
                                broadwidth, nbins, resolution)
          !
        END DO
      END DO
    END IF
  END SUBROUTINE


    !###added in 4/22/12
    SUBROUTINE allocate_spec(di_spectrum,sc_spectrum,sv_spectrum,ispin,cutoff,broadwidth,resolution)
        real(dp),dimension(:,:,:,:,:), allocatable                              :: di_spectrum
        real(dp),dimension(:,:,:,:,:), allocatable                              :: sc_spectrum
        real(dp),dimension(:,:,:), allocatable                                  :: sv_spectrum
        real(dp),                                               intent(in)      :: cutoff, broadwidth, resolution
        integer,                                                intent(in)      :: ispin


        integer                                                                 :: nbins

        if (allocated(di_spectrum) .or. allocated(sc_spectrum) .or. allocated(sv_spectrum) ) return

        nbins=nint(cutoff/broadwidth*resolution)
        allocate(di_spectrum(2,nbins,3,3,(ceiling(ispin/2.0)+1)))
        allocate(sc_spectrum(nbins,3,3,3,(ceiling(ispin/2.0)+1)))
        allocate(sv_spectrum(nbins,3,(ceiling(ispin/2.0)+1)))

        di_spectrum=0.0
        sc_spectrum=0.0
        sv_spectrum=0.0

    end SUBROUTINE




    subroutine energy_grid_ev(egrid,cutoff,broadwidth,resolution)
        real(dp), dimension(:), allocatable,                    intent(out)     :: egrid
        real(dp),                                               intent(in)      :: broadwidth, cutoff
        real(dp),                                               intent(in)      :: resolution

        integer                                                                 :: ibin


        allocate(egrid(int(cutoff/broadwidth*resolution)))


        forall (ibin=1:size(egrid))  egrid(ibin)=(ibin-1)*broadwidth/resolution*hartree_ev

    end subroutine

  ! --------------------------------------------------------------------------
  SUBROUTINE collect_spectrum(spectrum)
  ! --------------------------------------------------------------------------
    USE mp_global, ONLY : mp_barrier, mp_reduce_sum, mp_comm_all
    !
    REAL(dp), INTENT(inout) :: spectrum(:)
    !
    CALL mp_barrier(mp_comm_all)
    CALL flush(6)
    !
    CALL mp_reduce_sum(spectrum, mp_comm_all)
    CALL flush(6)
  END SUBROUTINE collect_spectrum

    !###added in 4/22/12
    SUBROUTINE collect_spectrum_di(di_spectrum)
        real(dp), dimension(:,:,:,:,:),                                 intent(inout)   :: di_spectrum

        integer                                                                         :: st, i, j

        write(6,*)"Gathering di ..."
        do st=1,size(di_spectrum,5)
            do j=1,3
                do i=1,3
                    call collect_spectrum(di_spectrum(1,:,i,j,st))
                    call collect_spectrum(di_spectrum(2,:,i,j,st))
                enddo
            enddo
        end do

        write(6,*)"Done di ..."
        if (size(di_spectrum,5)==1) then
            di_spectrum(1,:,:,:,:)=di_spectrum(1,:,:,:,:)+1
        else
            di_spectrum(1,:,:,:,:)=di_spectrum(1,:,:,:,:)+0.5
        end if


    end SUBROUTINE

    !###added in 4/22/12
    SUBROUTINE collect_spectrum_sc(sc_spectrum)
        real(dp), dimension(:,:,:,:,:),                                 intent(inout)   :: sc_spectrum

        integer                                                                         :: i, j, q, st

        write(6,*)"Gathering sc ..."
        do st=1,size(sc_spectrum,5)
            do q=1,3
                do j=1,3
                    do i=1,3
                        call collect_spectrum(sc_spectrum(:,i,j,q,st))
                    enddo
                enddo
            end do
        end do
        write(6,*)"Done sc ..."
    end SUBROUTINE


    SUBROUTINE collect_spectrum_sg(sg_spectrum)
        real(dp), dimension(:,:,:,:,:,:),                               intent(inout)   :: sg_spectrum

        integer                                                                         :: i, j, q, st

        write(6,*)"Gathering sg ..."
        do st=1,size(sg_spectrum,6)
            do q=1,3
                do j=1,3
                    do i=1,3
                        call collect_spectrum(sg_spectrum(:,i,j,q,1,st))
                        call collect_spectrum(sg_spectrum(:,i,j,q,2,st))
                    enddo
                enddo
            end do
        end do
        write(6,*)"Done sg ..."
    end SUBROUTINE


    !###added in 4/22/12
    SUBROUTINE collect_spectrum_sc2(presc_spectrum)
        real(dp), dimension(:,:,:,:,:),                                 intent(inout)   :: presc_spectrum

        integer                                                                         :: i, j, q, st

        write(6,*)"Gathering presc ..."
        do st=1,size(presc_spectrum,5)
            do q=1,3
                do j=1,3
                    do i=1,3
                        call collect_spectrum(presc_spectrum(:,i,j,q,st))
                    enddo
                enddo
            end do
        end do
        write(6,*)"Done presc ..."
    end SUBROUTINE

    !###added in 4/22/12
    SUBROUTINE collect_spectrum_sv(sv_spectrum)
        real(dp), dimension(:,:,:),                                   intent(inout)   :: sv_spectrum

        integer                                                                       :: i, j, st

        write(6,*)"Gathering sv ..."
        do st=1,size(sv_spectrum,3)
            do j=1,3
                call collect_spectrum(sv_spectrum(:,j,st))
            enddo
        end do
        write(6,*)"Done sv ..."
    end SUBROUTINE




    pure function delta_broaden(y,freq,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: delta_broaden


        integer                                                                 :: nwidth
        integer                                                                 :: i, eindex

        nwidth=4*resolution
        delta_broaden=0d0

        eindex=NINT(freq/broadwidth*resolution)

        forall (i=max(eindex-nwidth,1):min(eindex+nwidth,nbins)) delta_broaden(i)=y*(.5d0*broadwidth)/(pi*(((eindex-i)/resolution)**2+.25d0)*broadwidth**2)

    end function

    pure function pp_broaden(y,freq,broadwidth,nbins,resolution)
        real(dp),                                               intent(in)      :: y, freq, broadwidth
        integer,                                                intent(in)      :: nbins
        real(dp),                                               intent(in)      :: resolution

        real(dp), dimension(nbins)                                              :: pp_broaden


        integer                                                                 :: i, eindex

        pp_broaden=0d0

        eindex=NINT(freq/broadwidth*resolution)

        forall (i=1:nbins) pp_broaden(i)=y*(((eindex-i)/resolution)*broadwidth)/( (((eindex-i)/resolution)**2+.25d0)*broadwidth**2 )

    end function

    pure FUNCTION ndelta_broaden_value(freq_x,freq_0,broadwidth,resolution)
        real(dp),                                               intent(in)      :: freq_x, freq_0, broadwidth
        real(dp),                                               intent(in)      :: resolution
        real(dp)                                                                :: ndelta_broaden_value

        integer                                                                 :: eindex_0, eindex_x

       !eindex_0=NINT(freq_0/broadwidth*resolution)
        eindex_x=NINT((freq_0-freq_x)/broadwidth*resolution)
        ndelta_broaden_value=(.5d0*broadwidth)/(pi*(((eindex_x)/resolution)**2+.25d0)*broadwidth**2)
    end FUNCTION

    pure FUNCTION npp_broaden_value(freq_x,freq_0,broadwidth,resolution)
        real(dp),                                               intent(in)      :: freq_x, freq_0, broadwidth
        real(dp),                                               intent(in)      :: resolution
        real(dp)                                                                :: npp_broaden_value

        integer                                                                 :: eindex_0, eindex_x

       !eindex_0=NINT(freq_0/broadwidth*resolution)
        eindex_x=NINT((freq_0-freq_x)/broadwidth*resolution)
        npp_broaden_value=(((eindex_x)/resolution)*broadwidth)/( (((eindex_x)/resolution)**2+.25d0)*broadwidth**2 )
    end FUNCTION



    pure real(dp) function n_r(di)
        real(dp), dimension(2),                                     intent(in)  :: di

        n_r=sqrt(sqrt(di(1)**2+di(2)**2)+di(1))/sqrt(2d0)

    end function

    pure real(dp) function alpha(di,eng)
        real(dp), dimension(2),                                     intent(in)  :: di
        real(dp),                                                   intent(in)  :: eng

        alpha=4d0*pi*eng/(h_ev*c_e)*sqrt(sqrt(di(1)**2+di(2)**2)-di(1))/sqrt(2d0)

    end function

    subroutine write_sv_spectrum(egrid,sv,sv2)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:),                                   intent(in)  :: sv
        real(dp), dimension(:,:), optional,                         intent(in)  :: sv2

        character(len=16), dimension(3)                                         :: sv_file

        integer                                                                 :: q,ieng

        sv_file(1)="sv.X.dat"
        sv_file(2)="sv.Y.dat"
        sv_file(3)="sv.Z.dat"

        write(6,*)"Printing sv spectra..."
        call flush(6)

        do q=1,3
            open(unit=20,file=sv_file(q),form='formatted')
                if (present(sv2)) then
                    write(20,'(4E22.10)') (egrid(ieng),sv(ieng,q),sv2(ieng,q),sv(ieng,q)+sv2(ieng,q),ieng=1,size(egrid))
                else
                    write(20,'(2E22.10)') (egrid(ieng),sv(ieng,q),ieng=1,size(egrid))
                end if
            close(20)
        enddo
        write(6,*)"Printing sv spectra finished"
    end subroutine

    !###modified in 4/22/12
    SUBROUTINE write_sv_spectrum_kl(egrid,sv)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:),                                 intent(in)  :: sv

        character(len=16), dimension(3)                                         :: sv_file

        integer                                                                 :: q,ieng

        sv_file(1)="sv.X.dat"
        sv_file(2)="sv.Y.dat"
        sv_file(3)="sv.Z.dat"

        write(6,*)"Printing sv spectra..."
        call flush(6)

        do q=1,3
            open(unit=20,file=sv_file(q),form='formatted')
                if (size(sv,3)==2) then
                    write(20,'(4E22.10)') (egrid(ieng),sv(ieng,q,1),sv(ieng,q,2),0.5*(sv(ieng,q,1)+sv(ieng,q,2)),ieng=1,size(egrid))
                else
                    write(20,'(2E22.10)') (egrid(ieng),sv(ieng,q,1),ieng=1,size(egrid))
                end if
            close(20)
        enddo
        write(6,*)"Printing sv spectra finished"
    end subroutine


    !###modified in 4/22/12
    subroutine write_di_spectrum_kl(egrid,di)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: di

        character(len=16), dimension(3,3)                                       :: di_file
        character(len=16), dimension(3,3)                                       :: n_file
        character(len=16), dimension(3,3)                                       :: alpha_file

        integer                                                                 :: i,j,ieng

        di_file(1,1)="di.xx.dat"
        di_file(2,1)="di.yx.dat"
        di_file(3,1)="di.zx.dat"
        di_file(2,2)="di.yy.dat"
        di_file(3,2)="di.zy.dat"
        di_file(3,3)="di.zz.dat"
        n_file(1,1)="n.xx.dat"
        n_file(2,1)="n.yx.dat"
        n_file(3,1)="n.zx.dat"
        n_file(2,2)="n.yy.dat"
        n_file(3,2)="n.zy.dat"
        n_file(3,3)="n.zz.dat"
        alpha_file(1,1)="alpha.xx.dat"
        alpha_file(2,1)="alpha.yx.dat"
        alpha_file(3,1)="alpha.zx.dat"
        alpha_file(2,2)="alpha.yy.dat"
        alpha_file(3,2)="alpha.zy.dat"
        alpha_file(3,3)="alpha.zz.dat"

        write(6,*)"Printing di spectra..."
        call flush(6)

        do i=1,3
            do j=i,3
                open(unit=20,file=di_file(j,i),form='formatted')
                    if (size(di,5)==2) then
                        write(20,'(4E22.10)') (egrid(ieng),di(1,ieng,i,j,1),di(1,ieng,i,j,2),di(1,ieng,i,j,1)+di(1,ieng,i,j,2),ieng=1,size(egrid))
                        write(20,*)
                        write(20,'(4E22.10)') (egrid(ieng),di(2,ieng,i,j,1),di(2,ieng,i,j,2),di(2,ieng,i,j,1)+di(2,ieng,i,j,2),ieng=1,size(egrid))
                    else
                        write(20,'(2E22.10)') (egrid(ieng),di(1,ieng,i,j,1),ieng=1,size(egrid))
                        write(20,*)
                        write(20,'(2E22.10)') (egrid(ieng),di(2,ieng,i,j,1),ieng=1,size(egrid))
                    end if
                close(20)

                open(unit=20,file=n_file(j,i),form='formatted')
                    if (size(di,5)==2) then
                        write(20,'(4E22.10)') (egrid(ieng),n_r(di(:,ieng,i,j,1)),n_r(di(:,ieng,i,j,2)),n_r(di(:,ieng,i,j,1)+di(:,ieng,i,j,2)),ieng=1,size(egrid))
                    else
                        write(20,'(2E22.10)') (egrid(ieng),n_r(di(:,ieng,i,j,1)),ieng=1,size(egrid))
                    end if
                close(20)
                open(unit=20,file=alpha_file(j,i),form='formatted')
                    if (size(di,5)==2) then
                         write(20,'(4E22.10)') (egrid(ieng),alpha(di(:,ieng,i,j,1),egrid(ieng)),alpha(di(:,ieng,i,j,2),egrid(ieng)),alpha(di(:,ieng,i,j,1)+di(:,ieng,i,j,2),egrid(ieng)),ieng=1,size(egrid))
                    else
                         write(20,'(2E22.10)') (egrid(ieng),alpha(di(:,ieng,i,j,1),egrid(ieng)),ieng=1,size(egrid))
                    end if
                close(20)
            enddo
        enddo
        write(6,*)"Printing di spectra finished"
        write(6,*)"Printing n spectra finished"
        write(6,*)"Printing alpha spectra finished"

    end subroutine


    !###modified in 4/22/12
    subroutine write_sc_spectrum_kl(egrid,sc,di)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: sc, di

        character(len=16), dimension(3,3,3)                                     :: sc_file

        integer                                                                 :: i,j,q,ieng

        sc_file(1,1,1)="sc.xxX.dat"
        sc_file(2,1,1)="sc.yxX.dat"
        sc_file(3,1,1)="sc.zxX.dat"
        sc_file(2,2,1)="sc.yyX.dat"
        sc_file(3,2,1)="sc.zyX.dat"
        sc_file(3,3,1)="sc.zzX.dat"

        sc_file(1,1,2)="sc.xxY.dat"
        sc_file(2,1,2)="sc.yxY.dat"
        sc_file(3,1,2)="sc.zxY.dat"
        sc_file(2,2,2)="sc.yyY.dat"
        sc_file(3,2,2)="sc.zyY.dat"
        sc_file(3,3,2)="sc.zzY.dat"

        sc_file(1,1,3)="sc.xxZ.dat"
        sc_file(2,1,3)="sc.yxZ.dat"
        sc_file(3,1,3)="sc.zxZ.dat"
        sc_file(2,2,3)="sc.yyZ.dat"
        sc_file(3,2,3)="sc.zyZ.dat"
        sc_file(3,3,3)="sc.zzZ.dat"

        write(6,*)"Printing sc spectra..."
        call flush(6)

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=sc_file(j,i,q),form='formatted')
                        if (size(sc,5)==2) then
                            write(20,'(5E22.10)') (egrid(ieng),sc(ieng,i,j,q,1),sc(ieng,i,j,q,2),sc(ieng,i,j,q,1)+sc(ieng,i,j,q,2),(sc(ieng,i,j,q,1)+sc(ieng,i,j,q,2))*n_r(di(:,ieng,i,j,1)+di(:,ieng,i,j,2))/(2*c_e*eps_0*(di(1,ieng,i,j,1)+di(1,ieng,i,j,2))),ieng=1,size(egrid))
                        else
                            write(20,'(3E22.10)') (egrid(ieng),sc(ieng,i,j,q,1),sc(ieng,i,j,q,1)*n_r(di(:,ieng,i,j,1))/(2*c_e*eps_0*di(1,ieng,i,j,1)),ieng=1,size(egrid))
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing sc spectra finished"

    end subroutine

    SUBROUTINE write_sg_spectrum_kl(egrid,sg)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:,:),                           intent(in)  :: sg

        character(len=16), dimension(3,3,3)                                     :: sg_file

        integer                                                                 :: i,j,q,ieng

        sg_file(1,1,1)="sg.xxX.dat"
        sg_file(2,1,1)="sg.yxX.dat"
        sg_file(3,1,1)="sg.zxX.dat"
        sg_file(1,2,1)="sg.xyX.dat"
        sg_file(2,2,1)="sg.yyX.dat"
        sg_file(3,2,1)="sg.zyX.dat"
        sg_file(1,3,1)="sg.xzX.dat"
        sg_file(2,3,1)="sg.yzX.dat"
        sg_file(3,3,1)="sg.zzX.dat"

        sg_file(1,1,2)="sg.xxY.dat"
        sg_file(2,1,2)="sg.yxY.dat"
        sg_file(3,1,2)="sg.zxY.dat"
        sg_file(1,2,2)="sg.yyY.dat"
        sg_file(2,2,2)="sg.yyY.dat"
        sg_file(3,2,2)="sg.zyY.dat"
        sg_file(1,3,2)="sg.xzY.dat"
        sg_file(2,3,2)="sg.yzY.dat"
        sg_file(3,3,2)="sg.zzY.dat"

        sg_file(1,1,3)="sg.xxZ.dat"
        sg_file(2,1,3)="sg.yxZ.dat"
        sg_file(3,1,3)="sg.zxZ.dat"
        sg_file(1,2,3)="sg.xyZ.dat"
        sg_file(2,2,3)="sg.yyZ.dat"
        sg_file(3,2,3)="sg.zyZ.dat"
        sg_file(1,3,3)="sg.xzZ.dat"
        sg_file(2,3,3)="sg.yzZ.dat"
        sg_file(3,3,3)="sg.zzZ.dat"

        write(6,*)"Printing sg spectra..."
        call flush(6)

        do q=1,3
            do i=1,3
                do j=1,3
                    open(unit=20,file=sg_file(j,i,q),form='formatted')
                        if (size(sg,6)==2) then
                            write(20,'(4E22.10)') (egrid(ieng),sg(ieng,i,j,q,1,1),sg(ieng,i,j,q,1,2),sg(ieng,i,j,q,1,1)+sg(ieng,i,j,q,1,2),ieng=1,size(egrid))
                            write(20,*)
                            write(20,'(4E22.10)') (egrid(ieng),sg(ieng,i,j,q,2,1),sg(ieng,i,j,q,2,2),sg(ieng,i,j,q,2,1)+sg(ieng,i,j,q,2,2),ieng=1,size(egrid))
                        else
                            write(20,'(2E22.10)') (egrid(ieng),sg(ieng,i,j,q,1,1),ieng=1,size(egrid))
                            write(20,*)
                            write(20,'(2E22.10)') (egrid(ieng),sg(ieng,i,j,q,2,1),ieng=1,size(egrid))
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing sg spectra finished"

    end SUBROUTINE




    subroutine write_glass_spectrum(egrid,sc,di,sc2,di2)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:),                               intent(in)  :: sc,di
        real(dp), dimension(:,:,:,:), optional,                     intent(in)  :: sc2,di2

        character(len=16), dimension(3,3,3)                                     :: glass_file

        integer                                                                 :: i,j,q,ieng

        glass_file(1,1,1)="glass.xxX.dat"
        glass_file(2,1,1)="glass.yxX.dat"
        glass_file(3,1,1)="glass.zxX.dat"
        glass_file(2,2,1)="glass.yyX.dat"
        glass_file(3,2,1)="glass.zyX.dat"
        glass_file(3,3,1)="glass.zzX.dat"

        glass_file(1,1,2)="glass.xxY.dat"
        glass_file(2,1,2)="glass.yxY.dat"
        glass_file(3,1,2)="glass.zxY.dat"
        glass_file(2,2,2)="glass.yyY.dat"
        glass_file(3,2,2)="glass.zyY.dat"
        glass_file(3,3,2)="glass.zzY.dat"

        glass_file(1,1,3)="glass.xxZ.dat"
        glass_file(2,1,3)="glass.yxZ.dat"
        glass_file(3,1,3)="glass.zxZ.dat"
        glass_file(2,2,3)="glass.yyZ.dat"
        glass_file(3,2,3)="glass.zyZ.dat"
        glass_file(3,3,3)="glass.zzZ.dat"

        write(6,*)"Printing sc spectra..."
        call flush(6)

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=glass_file(j,i,q),form='formatted')
                        if (present(sc2) .and. present(di2)) then
                            write(20,'(4E22.10)') (egrid(ieng),sc(ieng,i,j,q)*n_r(di(:,ieng,i,j))/(2d0*di(1,ieng,i,j)*eps_0*c_e)/(alpha(di(:,ieng,i,j),egrid(ieng))+1.5d6), sc2(ieng,i,j,q)*n_r(di2(:,ieng,i,j))/(2d0*di2(1,ieng,i,j)*eps_0*c_e)/(alpha(di2(:,ieng,i,j),egrid(ieng))+1.5d6), (sc(ieng,i,j,q)+sc2(ieng,i,j,q))*n_r(di(:,ieng,i,j)+di2(:,ieng,i,j))/(2d0*(di(1,ieng,i,j)+di2(1,ieng,i,j))*eps_0*c_e)/(alpha(di(:,ieng,i,j)+di2(:,ieng,i,j),egrid(ieng))+1.5d6),ieng=1,size(egrid))
                        else
                            write(20,'(2E22.10)') (egrid(ieng),sc(ieng,i,j,q)*n_r(di(:,ieng,i,j))/(2d0*di(1,ieng,i,j)*eps_0*c_e)/(alpha(di(:,ieng,i,j),egrid(ieng))+1.5d6),ieng=1,size(egrid))
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing glass spectra finished"

    end subroutine



    !###modified in 4/22/12
    subroutine write_glass_spectrum_kl(egrid,sc,di)
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:),                               intent(in)  :: sc,di

        character(len=16), dimension(3,3,3)                                     :: glass_file

        integer                                                                 :: i,j,q,ieng
        real(dp)                                                                :: thick, absorp, scaler

        glass_file(1,1,1)="glass.xxX.dat"
        glass_file(2,1,1)="glass.yxX.dat"
        glass_file(3,1,1)="glass.zxX.dat"
        glass_file(2,2,1)="glass.yyX.dat"
        glass_file(3,2,1)="glass.zyX.dat"
        glass_file(3,3,1)="glass.zzX.dat"

        glass_file(1,1,2)="glass.xxY.dat"
        glass_file(2,1,2)="glass.yxY.dat"
        glass_file(3,1,2)="glass.zxY.dat"
        glass_file(2,2,2)="glass.yyY.dat"
        glass_file(3,2,2)="glass.zyY.dat"
        glass_file(3,3,2)="glass.zzY.dat"

        glass_file(1,1,3)="glass.xxZ.dat"
        glass_file(2,1,3)="glass.yxZ.dat"
        glass_file(3,1,3)="glass.zxZ.dat"
        glass_file(2,2,3)="glass.yyZ.dat"
        glass_file(3,2,3)="glass.zyZ.dat"
        glass_file(3,3,3)="glass.zzZ.dat"

        write(6,*)"Printing glass spectra..."
        call flush(6)

        thick=130d-9

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=glass_file(j,i,q),form='formatted')
                        if (size(di,5)==2) then
                            do ieng=1,size(egrid)
                                absorp=alpha(di(:,ieng,i,j,1)+di(:,ieng,i,j,2),egrid(ieng))
                                scaler=1-exp(-absorp*thick)
!                               if (absorp>=1e5) then
!                                   scaler=1.0
!                               else
!                                   scaler=1-exp(-absorp*thick)
!                               end if
                                write(20,'(4E22.10)') egrid(ieng),(sc(ieng,i,j,q,1)+sc(ieng,i,j,q,2))*n_r(di(:,ieng,i,j,1)+di(:,ieng,i,j,2))/(2d0*(di(1,ieng,i,j,1)+di(1,ieng,i,j,2))*eps_0*c_e)/(absorp+1.0)*scaler
                            end do
                        else
                            do ieng=1,size(egrid)
                                absorp=alpha(di(:,ieng,i,j,1),egrid(ieng))
                                scaler=1-exp(-absorp*thick)
!                               if (absorp>=1e5) then
!                                   scaler=1.0
!                               else
!                                   scaler=1-exp(-absorp*thick)
!                               end if
                                write(20,'(2E22.10)') egrid(ieng),sc(ieng,i,j,q,1)*n_r(di(:,ieng,i,j,1))/(2d0*di(1,ieng,i,j,1)*eps_0*c_e)/(absorp+1.0)*scaler
                            end do
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing glass spectra finished"

    end subroutine




!///Modify Glass Coefficient
    !###modified in 4/22/12
    SUBROUTINE calculate_glass_spectrum(egrid,pre_sc,di,cutoff,broadwidth,resolution)
        real(dp),                                                   intent(in)  :: broadwidth, cutoff, resolution
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: pre_sc
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: di

        character(len=20), dimension(3,3,3)                                     :: glass_file
        real(dp),dimension(nint(cutoff/broadwidth*resolution),3,3,3)            :: glass_spec, glass_spec2
        integer                                                                 :: i,j,q,ibin, ieng, nbins

        glass_file(1,1,1)="Modi_glass.xxX.dat"
        glass_file(2,1,1)="Modi_glass.yxX.dat"
        glass_file(3,1,1)="Modi_glass.zxX.dat"
        glass_file(2,2,1)="Modi_glass.yyX.dat"
        glass_file(3,2,1)="Modi_glass.zyX.dat"
        glass_file(3,3,1)="Modi_glass.zzX.dat"

        glass_file(1,1,2)="Modi_glass.xxY.dat"
        glass_file(2,1,2)="Modi_glass.yxY.dat"
        glass_file(3,1,2)="Modi_glass.zxY.dat"
        glass_file(2,2,2)="Modi_glass.yyY.dat"
        glass_file(3,2,2)="Modi_glass.zyY.dat"
        glass_file(3,3,2)="Modi_glass.zzY.dat"

        glass_file(1,1,3)="Modi_glass.xxZ.dat"
        glass_file(2,1,3)="Modi_glass.yxZ.dat"
        glass_file(3,1,3)="Modi_glass.zxZ.dat"
        glass_file(2,2,3)="Modi_glass.yyZ.dat"
        glass_file(3,2,3)="Modi_glass.zyZ.dat"
        glass_file(3,3,3)="Modi_glass.zzZ.dat"

        glass_spec=0
        glass_spec2=0
        nbins=nint(cutoff/broadwidth*resolution)

        do q=1,3
           do i=1,3
              do j=1,3
                 do ibin=1,nbins
                    if (size(pre_sc,5)==2 .and. size(di,5)==2) then
                       glass_spec(:,i,j,q)=glass_spec(:,i,j,q)+ delta_broaden_index( (pre_sc(ibin,i,j,q,1)+pre_sc(ibin,i,j,q,2))*n_r(di(:,ibin,i,j,1)+di(:,ibin,i,j,2))/(2d0*(di(1,ibin,i,j,1)+di(1,ibin,i,j,2))*eps_0*c_e)/(alpha(di(:,ibin,i,j,1)+di(:,ibin,i,j,2),egrid(ibin))+1.0),ibin,broadwidth,nbins,resolution )
                    else
                       glass_spec(:,i,j,q)=glass_spec(:,i,j,q)+ delta_broaden_index(  pre_sc(ibin,i,j,q,1)*n_r(di(:,ibin,i,j,1))/(2d0*di(1,ibin,i,j,1)*eps_0*c_e)/(alpha(di(:,ibin,i,j,1),egrid(ibin))+1), ibin,broadwidth,nbins,resolution)
                    end if
                 end do
              end do
           end do
        end do

        write(6,*)"Printing Modified glass spectra..."
        call flush(6)

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=glass_file(j,i,q),form='formatted')
                        if (size(pre_sc,5)==2 .and. size(di,5)==2) then
                            write(20,'(2E22.10)') (egrid(ieng),glass_spec(ieng,i,j,q),ieng=1,size(egrid))
                        else
                            write(20,'(2E22.10)') (egrid(ieng),glass_spec(ieng,i,j,q),ieng=1,size(egrid))
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing Modified glass spectra finished"

    end SUBROUTINE
!///Modify Glass Coefficient
    !###modified in 4/22/12
    SUBROUTINE calculate_glass_spectrum2(egrid,pre_sc,di,cutoff,broadwidth,resolution)
        real(dp),                                                   intent(in)  :: broadwidth, cutoff, resolution
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: pre_sc
        real(dp), dimension(:,:,:,:,:),                             intent(in)  :: di

        character(len=20), dimension(3,3,3)                                     :: glass_file
        real(dp),dimension(nint(cutoff/broadwidth*resolution),3,3,3)            :: glass_spec, glass_spec2
        real(dp)                                                                :: thick, absorp, scaler
        integer                                                                 :: i,j,q,ibin, ieng, nbins

        glass_file(1,1,1)="Thin_glass.xxX.dat"
        glass_file(2,1,1)="Thin_glass.yxX.dat"
        glass_file(3,1,1)="Thin_glass.zxX.dat"
        glass_file(2,2,1)="Thin_glass.yyX.dat"
        glass_file(3,2,1)="Thin_glass.zyX.dat"
        glass_file(3,3,1)="Thin_glass.zzX.dat"

        glass_file(1,1,2)="Thin_glass.xxY.dat"
        glass_file(2,1,2)="Thin_glass.yxY.dat"
        glass_file(3,1,2)="Thin_glass.zxY.dat"
        glass_file(2,2,2)="Thin_glass.yyY.dat"
        glass_file(3,2,2)="Thin_glass.zyY.dat"
        glass_file(3,3,2)="Thin_glass.zzY.dat"

        glass_file(1,1,3)="Thin_glass.xxZ.dat"
        glass_file(2,1,3)="Thin_glass.yxZ.dat"
        glass_file(3,1,3)="Thin_glass.zxZ.dat"
        glass_file(2,2,3)="Thin_glass.yyZ.dat"
        glass_file(3,2,3)="Thin_glass.zyZ.dat"
        glass_file(3,3,3)="Thin_glass.zzZ.dat"

        glass_spec=0
        glass_spec2=0
        nbins=nint(cutoff/broadwidth*resolution)

        ! unit of nm
        thick=130.0e-9

        do q=1,3
           do i=1,3
              do j=1,3
                 do ibin=1,nbins
                    if (size(pre_sc,5)==2 .and. size(di,5)==2) then
                       absorp=alpha(di(:,ibin,i,j,1)+di(:,ibin,i,j,2),egrid(ibin))
                       scaler=1-exp(-absorp*thick)
                       glass_spec(:,i,j,q) = glass_spec(:,i,j,q) &
                         + delta_broaden_index((pre_sc(ibin,i,j,q,1)+pre_sc(ibin,i,j,q,2)) &
                             * n_r(di(:,ibin,i,j,1)+di(:,ibin,i,j,2)) &
                             / (2d0*(di(1,ibin,i,j,1)+di(1,ibin,i,j,2))*eps_0*c_e)/(absorp+1.0)*scaler, &
                             ibin, broadwidth, nbins, resolution)
                    else
                       absorp=alpha(di(:,ibin,i,j,1),egrid(ibin))
                       scaler=1-exp(-absorp*thick)
                       glass_spec(:,i,j,q) = glass_spec(:,i,j,q) &
                           + delta_broaden_index(pre_sc(ibin,i,j,q,1)*n_r(di(:,ibin,i,j,1)) &
                               / (2d0*di(1,ibin,i,j,1)*eps_0*c_e)/(absorp+1.0)*scaler, &
                            ibin, broadwidth, nbins, resolution )
                    end if
                 end do
              end do
           end do
        end do

        write(6,*)"Printing Modified glass spectra..."
        call flush(6)

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=glass_file(j,i,q),form='formatted')
                        if (size(pre_sc,5)==2 .and. size(di,5)==2) then
                            write(20,'(2E20.10)') (egrid(ieng),glass_spec(ieng,i,j,q),ieng=1,size(egrid))
                        else
                            write(20,'(2E20.10)') (egrid(ieng),glass_spec(ieng,i,j,q),ieng=1,size(egrid))
                        end if
                    close(20)
                enddo
            enddo
        enddo
        write(6,*)"Printing Modified glass spectra finished"

    end SUBROUTINE



    SUBROUTINE spin_collect(ispin,sc_spec,pre_sc_spec,di_spec,sv_spec,sc_spec2,pre_sc_spec2,di_spec2,sv_spec2)
        USE mp_global, ONLY : proc_id

        integer,                                                   intent(in)     :: ispin
        real(dp), dimension(:,:), allocatable,                     intent(inout)  :: sv_spec
        real(dp), dimension(:,:), allocatable,                     intent(out)    :: sv_spec2
        real(dp), dimension(:,:,:,:), allocatable,                 intent(inout)  :: sc_spec
        real(dp), dimension(:,:,:,:), allocatable,                 intent(out)    :: sc_spec2
        real(dp), dimension(:,:,:,:), allocatable,                 intent(inout)  :: pre_sc_spec
        real(dp), dimension(:,:,:,:), allocatable,                 intent(out)    :: pre_sc_spec2
        real(dp), dimension(:,:,:,:), allocatable,                 intent(inout)  :: di_spec
        real(dp), dimension(:,:,:,:), allocatable,                 intent(out)    :: di_spec2

        if ( ispin == 2 ) return

        if ( ispin == 1 ) then


                allocate(di_spec2(size(di_spec,1),size(di_spec,2),size(di_spec,3),size(di_spec,4)))
                allocate(sc_spec2(size(sc_spec,1),size(sc_spec,2),size(sc_spec,3),size(sc_spec,4)))
                allocate(pre_sc_spec2(size(pre_sc_spec,1),size(pre_sc_spec,2),size(pre_sc_spec,3),size(pre_sc_spec,4)))
                allocate(sv_spec2(size(sv_spec,1),size(sv_spec,2)))

                sc_spec2=sc_spec
                pre_sc_spec2=pre_sc_spec
                di_spec2=di_spec
                sv_spec2=sv_spec

!               write(6,*) proc_id, sc_spec2(50,1,1,1), ispin



            if ( allocated(sc_spec) ) deallocate(sc_spec)
            if ( allocated(pre_sc_spec) ) deallocate(pre_sc_spec)
            if ( allocated(di_spec) ) deallocate(di_spec)
            if ( allocated(sv_spec) ) deallocate(sv_spec)

        end if

    if (proc_id==0) write(6,*) "!!!!!!!!!spectrum matrices are copied and spectrum has been finished"

    END SUBROUTINE




    subroutine write_2sc_spectrum(egrid,sc,sc2)
        USE mp_global, ONLY : proc_id
        real(dp), dimension(:),                                     intent(in)  :: egrid
        real(dp), dimension(:,:,:,:),                               intent(in)  :: sc
        real(dp), dimension(:,:,:,:),                               intent(in)  :: sc2

        character(len=16), dimension(3,3,3)                                     :: sc_file

        integer                                                                 :: i,j,q,ieng

        sc_file(1,1,1)="sc.xxX.dat"
        sc_file(2,1,1)="sc.yxX.dat"
        sc_file(3,1,1)="sc.zxX.dat"
        sc_file(2,2,1)="sc.yyX.dat"
        sc_file(3,2,1)="sc.zyX.dat"
        sc_file(3,3,1)="sc.zzX.dat"

        sc_file(1,1,2)="sc.xxY.dat"
        sc_file(2,1,2)="sc.yxY.dat"
        sc_file(3,1,2)="sc.zxY.dat"
        sc_file(2,2,2)="sc.yyY.dat"
        sc_file(3,2,2)="sc.zyY.dat"
        sc_file(3,3,2)="sc.zzY.dat"

        sc_file(1,1,3)="sc.xxZ.dat"
        sc_file(2,1,3)="sc.yxZ.dat"
        sc_file(3,1,3)="sc.zxZ.dat"
        sc_file(2,2,3)="sc.yyZ.dat"
        sc_file(3,2,3)="sc.zyZ.dat"
        sc_file(3,3,3)="sc.zzZ.dat"

        write(6,*)"Printing sc spectra..."
        write(6,*) proc_id, sc(50,1,1,1), sc2(50,1,1,1)
        call flush(6)

        do q=1,3
            do i=1,3
                do j=i,3
                    open(unit=20,file=sc_file(j,i,q),form='formatted')
                        write(20,'(4E20.10)') (egrid(ieng),sc(ieng,i,j,q),sc2(ieng,i,j,q),sc(ieng,i,j,q)+sc2(ieng,i,j,q),ieng=1,size(egrid))
                    close(20)
                enddo
            enddo
        enddo

    end subroutine


end module
