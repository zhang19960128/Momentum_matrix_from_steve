module op_tools
    use flinal_tools
    use es_tools
    implicit none  

    type scalar_op
        complex(dpc), dimension(:,:), allocatable                           :: s
    end type

    type vector_op
        complex(dpc), dimension(:,:,:), allocatable                         :: v
    end type
    
    type vector_op_two
        complex(dpc), dimension(:,:), allocatable                           :: v
    end type

    type sc_op
        type (scalar_op), dimension(:), pointer                                 :: s_dk
        type (vector_op), dimension(:), pointer                                 :: p_dk
        type (vector_op), pointer                                               :: p
        real(dp)                                                            :: energy=0d0, filling=0d0
        integer, dimension(2,2)                                             :: band_range=0
    end type



    type kpt_ops
        type(sc_op), dimension(:,:),     allocatable                       :: block
        integer,                         dimension(-3:3)                   :: kdiff
        integer                                                            :: nblocks
        real(dp), dimension(3)                                             :: k, dk
    end type







interface dot_prod
    module procedure dot_prod_operator
end interface

interface transform
    module procedure transform_op_rvec, transform_op_rmat, transform_op_r3mat, transform_op_cvec
end interface

contains

!! operator functions!!!!!!!!!!!
  ! function op_blocks(es,ikpt,cutoff,filename)
  !     type(kpt_ops)                                                       :: op_blocks
  !     type(es_data),pointer                                               :: es        
  !     integer,                                       intent(in)           :: ikpt
  !     character(len=256),                            intent(in), optional :: filename

  !     
  !     type(wavefunction), dimension(:), allocatable                       :: dk_block_f,dk_block_b
  !     integer                                                             :: ibnd,iblock,jblock,q, i
  !     real(dp)                                                            :: cutoff
  !     
! !      if (present(filename)) then

! !      endif

  !     write(6,*) proc_id, "computing operator elements for kpt ... ", ikpt, es%kpt(ikpt)%nblock
  !     
  !     

  !     op_blocks%kdiff=es%kpt(ikpt)%kdiff
  !     op_blocks%k=es%kpt(ikpt)%g%k
  !     op_blocks%dk=1d0/es%kptdim

  !     op_blocks%nblocks=es%kpt(ikpt)%nblock
  !     allocate(op_blocks%block(es%kpt(ikpt)%nblock, es%kpt(ikpt)%nblock))

  !     do jblock=1, es%kpt(ikpt)%nblock
  !         do iblock=1, es%kpt(ikpt)%nblock
  !             op_blocks%block(iblock,jblock)%band_range(:,1)=es%kpt(ikpt)%block(iblock)%band_range
  !             op_blocks%block(iblock,jblock)%band_range(:,2)=es%kpt(ikpt)%block(jblock)%band_range

  !             op_blocks%block(iblock,jblock)%energy=es%kpt(ikpt)%block(jblock)%energy-es%kpt(ikpt)%block(iblock)%energy
  !             op_blocks%block(iblock,jblock)%filling=es%kpt(ikpt)%block(jblock)%filling-es%kpt(ikpt)%block(iblock)%filling
  !         enddo
  !     enddo
  !       

  !     do jblock=1, es%kpt(ikpt)%nblock
  !         allocate(op_blocks%block(jblock,jblock)%s_dk(-3:3))
  !         do iblock=1, es%kpt(ikpt)%nblock
  !             if (abs(op_blocks%block(iblock,jblock)%energy)<cutoff .and. abs(op_blocks%block(iblock,jblock)%filling) > 0d0) then
  !                 allocate(op_blocks%block(iblock,jblock)%p)
  !                 allocate(op_blocks%block(iblock,jblock)%p_dk(-3:3))
  !                 op_blocks%block(iblock,jblock)%p=p_block(es%kpt(ikpt)%block(iblock)%wf,es%kpt(ikpt)%block(jblock)%wf,es%lattice%gprimd)
! !                 write(6,*) proc_id, ikpt, op_blocks%block(iblock,jblock)%p%v(1,1,1)
  !             else
  !                 nullify(op_blocks%block(iblock,jblock)%p,op_blocks%block(iblock,jblock)%p_dk)
  !             endif
  !         enddo
  !     enddo

  !     do jblock=1, es%kpt(ikpt)%nblock
  !         do q=1,3
  !             call get_dk_block(dk_block_f,es%kpt(ikpt)%block(jblock),es%kpt(es%kpt(ikpt)%kdiff(q))%block)
  !             call get_dk_block(dk_block_b,es%kpt(ikpt)%block(jblock),es%kpt(es%kpt(ikpt)%kdiff(-q))%block)
  ! 
  !             op_blocks%block(jblock,jblock)%s_dk(q)=s_block(es%kpt(ikpt)%block(jblock)%wf,dk_block_f)
  !             op_blocks%block(jblock,jblock)%s_dk(-q)=s_block(es%kpt(ikpt)%block(jblock)%wf,dk_block_b)

  !             do iblock=1, es%kpt(ikpt)%nblock
  !                 if (abs(op_blocks%block(iblock,jblock)%energy)<cutoff .and. abs(op_blocks%block(iblock,jblock)%filling) > 0d0) then
  !                     op_blocks%block(iblock,jblock)%p_dk(q)=p_block(es%kpt(ikpt)%block(iblock)%wf,dk_block_f,es%lattice%gprimd)
  !                     op_blocks%block(iblock,jblock)%p_dk(-q)=p_block(es%kpt(ikpt)%block(iblock)%wf,dk_block_b,es%lattice%gprimd)
  !                 endif
  !             enddo      
  !         enddo
  !     enddo

  !     if (allocated(dk_block_f)) then
  !         do ibnd=lbound(dk_block_f,1),ubound(dk_block_f,1)
  !             deallocate(dk_block_f(ibnd)%g,dk_block_f(ibnd)%coeff)
  !         enddo
  !         deallocate(dk_block_f)
  !     endif
  !     if (allocated(dk_block_b)) then
  !         do ibnd=lbound(dk_block_b,1),ubound(dk_block_b,1)
  !             deallocate(dk_block_b(ibnd)%g,dk_block_b(ibnd)%coeff)
  !         enddo
  !         deallocate(dk_block_b)
  !     endif

  ! end function

    !###modified in 4/22/12 
    FUNCTION op_blocks_kl(es_kl,cutoff,sg,filename)
        USE mp_global, ONLY : proc_id
        type(kpt_ops)                                                       :: op_blocks_kl
        type(es_header),pointer                                             :: es_kl        
        logical,                                       intent(in)           :: sg
        character(len=256),                            intent(in), optional :: filename

        
        type(wavefunction), dimension(:), allocatable                       :: dk_block_f,dk_block_b
        integer                                                             :: ibnd,iblock,jblock,q, i, vband
        real(dp)                                                            :: cutoff
        

        if (.not. allocated(es_kl%kpt(0)%block) ) return


        if (proc_id==0) write(6,*) "computing operator elements for kpt ... " !

        op_blocks_kl%kdiff=es_kl%kpt(0)%kdiff
        op_blocks_kl%k=es_kl%kpt(0)%g%k
        op_blocks_kl%dk=1d0/es_kl%kptdim

        op_blocks_kl%nblocks=es_kl%kpt(0)%nblock
        allocate(op_blocks_kl%block(es_kl%kpt(0)%nblock, es_kl%kpt(0)%nblock))

        do jblock=1, es_kl%kpt(0)%nblock
            do iblock=1, es_kl%kpt(0)%nblock
                op_blocks_kl%block(iblock,jblock)%band_range(:,1)=es_kl%kpt(0)%block(iblock)%band_range
                op_blocks_kl%block(iblock,jblock)%band_range(:,2)=es_kl%kpt(0)%block(jblock)%band_range

                op_blocks_kl%block(iblock,jblock)%energy=es_kl%kpt(0)%block(jblock)%energy-es_kl%kpt(0)%block(iblock)%energy
                op_blocks_kl%block(iblock,jblock)%filling=es_kl%kpt(0)%block(jblock)%filling-es_kl%kpt(0)%block(iblock)%filling
                if (abs(op_blocks_kl%block(iblock,jblock)%filling)<=eq_tol) op_blocks_kl%block(iblock,jblock)%filling=0d0
            enddo
        enddo
          
        vband=es_kl%kpt(0)%vband
        do jblock=1, es_kl%kpt(0)%nblock
            allocate(op_blocks_kl%block(jblock,jblock)%s_dk(-3:3))
            do iblock=1, es_kl%kpt(0)%nblock
    !           if (abs(op_blocks_kl%block(iblock,jblock)%energy)<cutoff .and. abs(op_blocks_kl%block(iblock,jblock)%filling) > 0d0) then
     !          allocate(op_blocks_kl%block(iblock,jblock)%p)
     !          op_blocks_kl%block(iblock,jblock)%p=p_block(es_kl%kpt(0)%block(iblock)%wf,es_kl%kpt(0)%block(jblock)%wf,es_kl%lattice%gprimd)
    !           else
                if (.not. sg) then
                    IF ( abs(op_blocks_kl%block(iblock,jblock)%energy) < cutoff &
                        .and. abs(op_blocks_kl%block(iblock,jblock)%filling) > 0d0 ) then
                        allocate(op_blocks_kl%block(iblock,jblock)%p)
                        op_blocks_kl%block(iblock,jblock)%p &
                            = p_block(es_kl%kpt(0)%block(iblock)%wf, &
                                      es_kl%kpt(0)%block(jblock)%wf, &
                                      es_kl%lattice%gprimd)
                    else
                        nullify(op_blocks_kl%block(iblock,jblock)%p_dk)
                    endif
                else
                    allocate(op_blocks_kl%block(iblock,jblock)%p)
                    op_blocks_kl%block(iblock,jblock)%p = &
                        p_block(es_kl%kpt(0)%block(iblock)%wf,es_kl%kpt(0)%block(jblock)%wf,es_kl%lattice%gprimd)
                end if

                if ( abs(op_blocks_kl%block(iblock,jblock)%energy) <= cutoff &
                   .and. abs(op_blocks_kl%block(iblock,jblock)%filling) > 0d0 ) then
                    allocate(op_blocks_kl%block(iblock,jblock)%p_dk(-3:3))
                else
                    nullify(op_blocks_kl%block(iblock,jblock)%p_dk)
                endif
    !           end if
            enddo
        enddo

        if (es_kl%kpts_type=='path') return

        do jblock=1, es_kl%kpt(0)%nblock
            do q=1,3
                call get_dk_block(dk_block_f,es_kl%kpt(0)%block(jblock),es_kl%kpt(q)%block)
                call get_dk_block(dk_block_b,es_kl%kpt(0)%block(jblock),es_kl%kpt(-q)%block)
    
                op_blocks_kl%block(jblock,jblock)%s_dk(q)=s_block(es_kl%kpt(0)%block(jblock)%wf,dk_block_f)
                op_blocks_kl%block(jblock,jblock)%s_dk(-q)=s_block(es_kl%kpt(0)%block(jblock)%wf,dk_block_b)

                do iblock=1, es_kl%kpt(0)%nblock
                    if ( abs(op_blocks_kl%block(iblock,jblock)%energy) <= cutoff &
                       .and. abs(op_blocks_kl%block(iblock,jblock)%filling) > 0d0) then
                        op_blocks_kl%block(iblock,jblock)%p_dk(q) = &
                            p_block(es_kl%kpt(0)%block(iblock)%wf, dk_block_f, es_kl%lattice%gprimd)
                        op_blocks_kl%block(iblock,jblock)%p_dk(-q) = &
                            p_block(es_kl%kpt(0)%block(iblock)%wf, dk_block_b, es_kl%lattice%gprimd)
                    endif
                enddo      
            enddo
        enddo

        if (allocated(dk_block_f)) then
            do ibnd=lbound(dk_block_f,1),ubound(dk_block_f,1)
                deallocate(dk_block_f(ibnd)%g,dk_block_f(ibnd)%coeff)
            enddo
            deallocate(dk_block_f)
        endif
        if (allocated(dk_block_b)) then
            do ibnd=lbound(dk_block_b,1),ubound(dk_block_b,1)
                deallocate(dk_block_b(ibnd)%g,dk_block_b(ibnd)%coeff)
            enddo
            deallocate(dk_block_b)
        endif

    end function


    function p_block(wf1,wf2,gprimd)
    !sets transition dipole matrix elements
        type(vector_op)                                                     :: p_block

        type (wavefunction), dimension(:),               intent(in)         :: wf1,wf2
        real(dp), dimension(3,3)                                            :: gprimd

    
        integer                                                             :: ibnd,jbnd

        allocate(p_block%v(size(wf1),size(wf2),3))
    
        forall (ibnd=1:size(wf1),jbnd=1:size(wf2))    p_block%v(ibnd,jbnd,:)=p_overlap(wf1(ibnd),wf2(jbnd))
!        forall (ibnd=1:size(wf1),jbnd=1:size(wf2))    p_block%v(ibnd,jbnd,:)=matmul(gprimd,p_overlap(wf1(ibnd),wf2(jbnd)))
    end function

    function inner_block(wf1,wf2,gprimd)
    !sets transition dipole matrix elements
        type(vector_op_two)                                                     :: inner_block
        type (wavefunction), dimension(:),               intent(in)         :: wf1,wf2
        real(dp), dimension(3,3)                                            :: gprimd

    
        integer                                                             :: ibnd,jbnd

        allocate(inner_block%v(size(wf1),size(wf2)))
    
        forall (ibnd=1:size(wf1),jbnd=1:size(wf2))    inner_block%v(ibnd,jbnd)=inner_prod(wf1(ibnd),wf2(jbnd))
!        forall (ibnd=1:size(wf1),jbnd=1:size(wf2))    p_block%v(ibnd,jbnd,:)=matmul(gprimd,p_overlap(wf1(ibnd),wf2(jbnd)))
    end function


    function s_block(wf1,wf2)
        type(scalar_op)                                                     :: s_block

        type (wavefunction), dimension(:),               intent(in)         :: wf1,wf2


        integer                                                             :: ibnd,jbnd

        allocate(s_block%s(size(wf1),size(wf2)))
    
        forall (ibnd=1:size(wf1),jbnd=1:size(wf2))    s_block%s(ibnd,jbnd)=inner_prod(wf1(ibnd),wf2(jbnd))
    end function

    

        



    

        


!    function compute_matrix(es,ikpt)
!        type(operators)                                                     :: compute_matrix
!        type(es_data)                                                       :: es        
!        integer,                                       intent(in)           :: ikpt
!
!        integer                                                             :: iblock,q
!
!        write(6,*) "computing operator elements for kpt ... ", ikpt
!        call flush(6)
!
!        allocate(compute_matrix%s_kdiff(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3), &
!                 compute_matrix%p_kdiff(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3,3), &
!                 compute_matrix%p(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3))
!        compute_matrix%p=momenta(es,ikpt)
!        compute_matrix%s_kdiff=overlap_grads(es,ikpt)
!        compute_matrix%p_kdiff=momentum_phase_grads(es,ikpt)
!
!
!        compute_matrix%kcoord=es%kpt(ikpt)%kcoord 
!        compute_matrix%dk=1d0/es%kptdim
!
!        compute_matrix%nblocks=size(es%kpt(ikpt)%block)
!        allocate(compute_matrix%block(size(es%kpt(ikpt)%block)))
!        compute_matrix%block=es%kpt(ikpt)%block
!
!        allocate(compute_matrix%dblock(2,size(es%kpt(ikpt)%block),-3:3))
!        do q=-3,3
!            do iblock=1,size(compute_matrix%block)
!                compute_matrix%dblock(:,iblock,q)=get_dk_range(compute_matrix%block(iblock),es%kpt(es%kpt(ikpt)%kdiff(q))%block)
!            enddo
!        enddo
!        
!    end function
!
!
!
!
!!!p_kdiff
!    function overlap_grads(es,ikpt)
!        type(es_data)                                                       :: es
!        integer,                                            intent(in)      :: ikpt
!
!        complex(dpc), dimension(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3)    :: overlap_grads
!
!        integer                                                             :: ibnd,jbnd,q
!        type(wavefunction)                                                  :: wf_temp
!        integer, dimension(3,3)                                             :: uvec
!
!        uvec=iidentity(3)
!
!        overlap_grads=0d0
!        do q=1,3        
!
!            ibnd=ubound(es%kpt(ikpt)%band,1)
!            do while (es%kpt(ikpt)%band(es%kpt(ikpt)%vband+1)%energy-es%kpt(ikpt)%band(ibnd)%energy<es%freq_cutoff)
!                if (allocated(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf%coeff)) then
!                    if (es%kpt(ikpt)%kcoord(q)-es%kpt(es%kpt(ikpt)%kdiff(q))%kcoord(q)>0d0) then
!                        wf_temp=reorder_basis(bz_wrap(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf,uvec(:,q)),es%kpt(ikpt)%g)
!                    else
!                        wf_temp=reorder_basis(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf,es%kpt(ikpt)%g)
!                    endif
!
!                    jbnd=ubound(es%kpt(ikpt)%band,1)           
!                    do while (es%kpt(ikpt)%band(es%kpt(ikpt)%vband+1)%energy-es%kpt(ikpt)%band(jbnd)%energy<es%freq_cutoff)
!                        overlap_grads(jbnd,ibnd,q)=inner_prod(es%kpt(ikpt)%band(jbnd)%wf,wf_temp)
!
!                        jbnd=jbnd-1
!                    enddo
!
!                    deallocate(wf_temp%coeff,wf_temp%g)
!                endif
!
!                ibnd=ibnd-1
!            enddo
!
!        enddo
!        
!        
!    end function
!   
!!!p
!    function momenta(es,ikpt)
!        type(es_data)                                                       :: es
!        integer,                                            intent(in)      :: ikpt
!
!        complex(dpc), dimension(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3)    :: momenta
!
!        integer                                                             :: ibnd,jbnd
!        
!        
!        momenta=0d0
!
!        ibnd=es%kpt(ikpt)%vband
!        do while (es%kpt(ikpt)%band(es%kpt(ikpt)%vband+1)%energy-es%kpt(ikpt)%band(ibnd)%energy<es%freq_cutoff)
!
!            jbnd=es%kpt(ikpt)%vband+1           
!            do while (es%kpt(ikpt)%band(jbnd)%energy-es%kpt(ikpt)%band(ibnd)%energy<es%freq_cutoff)
!                momenta(jbnd,ibnd,:)=matmul(es%lattice%gprimd,p_overlap(es%kpt(ikpt)%band(jbnd)%wf,es%kpt(ikpt)%band(ibnd)%wf))
!!                momenta(jbnd,ibnd,:)=p_overlap(es%kpt(ikpt)%band(jbnd)%wf,es%kpt(ikpt)%band(ibnd)%wf)
!                momenta(ibnd,jbnd,:)=conjg(momenta(jbnd,ibnd,:))
!
!                jbnd=jbnd+1
!            enddo
!
!            ibnd=ibnd-1
!        enddo
!
!    end function
!
!
!!!s_kdiff
!    function momentum_phase_grads(es,ikpt)
!        type(es_data)                                                       :: es
!        integer,                                            intent(in)      :: ikpt
!
!        complex(dpc), dimension(es%kpt(ikpt)%nband,es%kpt(ikpt)%nband,3,3)  :: momentum_phase_grads
!
!        integer                                                             :: ibnd,jbnd,q
!        type(wavefunction)                                                  :: wf_temp
!        integer, dimension(3,3)                                             :: uvec
!
!        uvec=iidentity(3)
!        
!        momentum_phase_grads=0d0
!        do q=1,3
!
!            ibnd=ubound(es%kpt(ikpt)%band,1)
!            do while (es%kpt(ikpt)%band(es%kpt(ikpt)%vband+1)%energy-es%kpt(ikpt)%band(ibnd)%energy<es%freq_cutoff)
!                if (allocated(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf%coeff)) then
!                    if (es%kpt(ikpt)%kcoord(q)-es%kpt(es%kpt(ikpt)%kdiff(q))%kcoord(q)>0d0) then
!                        wf_temp=reorder_basis(bz_wrap(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf,uvec(:,q)),es%kpt(ikpt)%g)
!                    else
!                        wf_temp=reorder_basis(es%kpt(es%kpt(ikpt)%kdiff(q))%band(ibnd)%wf,es%kpt(ikpt)%g)
!                    endif
!
!                jbnd=ubound(es%kpt(ikpt)%band,1)
!                    do while (es%kpt(ikpt)%band(es%kpt(ikpt)%vband+1)%energy-es%kpt(ikpt)%band(jbnd)%energy<es%freq_cutoff)
!                        momentum_phase_grads(jbnd,ibnd,:,q)=matmul(es%lattice%gprimd,p_overlap(es%kpt(ikpt)%band(jbnd)%wf,wf_temp))
!!                        momentum_phase_grads(jbnd,ibnd,:,q)=p_overlap(es%kpt(ikpt)%band(jbnd)%wf,wf_temp)
!
!                        jbnd=jbnd-1
!                    enddo
!
!                    deallocate(wf_temp%coeff,wf_temp%g)
!                endif
!
!                ibnd=ibnd-1
!            enddo
!        enddo
!
!    end function


    subroutine get_dk_block(wfblock,nblock, kblocks)
        type(wavefunction), dimension(:), allocatable                           :: wfblock

        type(block),                                            intent(in)      :: nblock
        type(block), dimension(:),                              intent(in)      :: kblocks



        type(wavefunction)                                                      :: bzw_temp
        integer                                                                 :: i,j,iblock,jblock, ibnd
        integer, dimension(3)                                                   :: kshift


        iblock=1
        do while(nblock%band_range(1) > kblocks(iblock)%band_range(2))
            iblock=iblock+1
        enddo

        jblock=iblock
        do while(nblock%band_range(2) > kblocks(jblock)%band_range(2))
            jblock=jblock+1
        enddo

        

        if (allocated(wfblock)) then
            do ibnd=lbound(wfblock,1),ubound(wfblock,1)
                deallocate(wfblock(ibnd)%g,wfblock(ibnd)%coeff)
            enddo
            deallocate(wfblock)
        endif

        allocate(wfblock(kblocks(iblock)%band_range(1):kblocks(jblock)%band_range(2)))

        if (any(abs(nblock%wf(ubound(nblock%wf,1))%g%k- kblocks(jblock)%wf(ubound(kblocks(jblock)%wf,1))%g%k) > .5)) then
            kshift=nint(nblock%wf(ubound(nblock%wf,1))%g%k- kblocks(jblock)%wf(ubound(kblocks(jblock)%wf,1))%g%k)
            do i=iblock,jblock
                do j=kblocks(i)%band_range(1),kblocks(i)%band_range(2)
                    bzw_temp=bz_wrap(kblocks(i)%wf(j) ,kshift)
                    wfblock(j)=reorder_basis(bzw_temp,nblock%wf(ubound(nblock%wf,1))%g)
                    deallocate(bzw_temp%g,bzw_temp%coeff)
                enddo
            enddo
        else
            do i=iblock,jblock
                do j=kblocks(i)%band_range(1),kblocks(i)%band_range(2)
                    wfblock(j)=reorder_basis(kblocks(i)%wf(j) ,nblock%wf(ubound(nblock%wf,1))%g)
                enddo
            enddo
        endif

    end subroutine

    function get_dk_range(nblock,kblocks)
        integer, dimension(2)                                                   :: get_dk_range

        type(block),                                            intent(in)      :: nblock
        type(block), dimension(:),                              intent(in)      :: kblocks
        

        integer                                                                 :: iblock        
                
      

        iblock=1
        do while(nblock%band_range(1) > kblocks(iblock)%band_range(2))
            iblock=iblock+1
        enddo
        get_dk_range(1)=kblocks(iblock)%band_range(1)

        do while(nblock%band_range(2) > kblocks(iblock)%band_range(2))
            iblock=iblock+1
        enddo
        get_dk_range(2)=kblocks(iblock)%band_range(2)

    end function


!!!!!!shift current functions!!!!!!!!!!!!!!

    function transform_op_rvec(op,lat_vec)
        real(dp), dimension(:,:,:),                             intent(in)      :: op
        real(dp), dimension(3,3),                               intent(in)      :: lat_vec

        real(dp), dimension(size(op,1),size(op,2),3)                            :: transform_op_rvec
        integer                                                                 :: i

        if (size(op,3) /= 3) then; write(6,*) "operator not rvec op"; stop; endif 
    
        forall(i=1:3)
            transform_op_rvec(:,:,i)=lat_vec(i,1)*op(:,:,1) + lat_vec(i,2)*op(:,:,2) + lat_vec(i,3)*op(:,:,3) 
        end forall
    end function


    function transform_op_rmat(op,lat_vec)
        real(dp), dimension(:,:,:,:),                           intent(in)      :: op
        real(dp), dimension(3,3),                               intent(in)      :: lat_vec

        real(dp), dimension(size(op,1),size(op,2),3,3)                          :: transform_op_rmat
        integer                                                                 :: i,j

        if (size(op,3) /= 3 .or. size(op,4) /= 3) then; write(6,*) "operator not rmat op"; stop; endif 
    
        forall(i=1:3,j=1:3)
            transform_op_rmat(:,:,i,j) &
                = lat_vec(i,1)*lat_vec(j,1)*op(:,:,1,1) &
                + lat_vec(i,2)*lat_vec(j,1)*op(:,:,2,1) &
                + lat_vec(i,3)*lat_vec(j,1)*op(:,:,3,1) &
                + lat_vec(i,1)*lat_vec(j,2)*op(:,:,1,2) &
                + lat_vec(i,2)*lat_vec(j,2)*op(:,:,2,2) &
                + lat_vec(i,3)*lat_vec(j,2)*op(:,:,3,2) &
                + lat_vec(i,1)*lat_vec(j,3)*op(:,:,1,3) &
                + lat_vec(i,2)*lat_vec(j,3)*op(:,:,2,3) &
                + lat_vec(i,3)*lat_vec(j,3)*op(:,:,3,3) 
        end forall
    end function



    FUNCTION transform_op_r3mat(op,lat_vec)
        real(dp), dimension(:,:,:,:,:),                         intent(in)      :: op
        real(dp), dimension(3,3),                               intent(in)      :: lat_vec

        real(dp), dimension(size(op,1),size(op,2),3,3,3)                        :: transform_op_r3mat
        integer                                                                 :: i,j

        if (size(op,3) /= 3 .or. size(op,4) /= 3 .or. size(op,5) /= 3) then; write(6,*) "operator not r3mat op"; stop; endif 

        do i=1,3
            transform_op_r3mat(:,:,:,:,i)=lat_vec(i,1)*transform_op_rmat(op(:,:,:,:,1),lat_vec) &
                                         +lat_vec(i,2)*transform_op_rmat(op(:,:,:,:,2),lat_vec) &
                                         +lat_vec(i,3)*transform_op_rmat(op(:,:,:,:,3),lat_vec)
        end do
    end FUNCTION
    
    function transform_op_cvec(op,lat_vec)
        complex(dpc), dimension(:,:,:),                         intent(in)      :: op
        real(dp), dimension(3,3),                               intent(in)      :: lat_vec

        complex(dpc), dimension(size(op,1),size(op,2),3)                        :: transform_op_cvec
        integer                                                                 :: i

        if (size(op,3) /= 3) then; write(6,*) "operator not cvec op"; stop; endif 
    
        forall(i=1:3)
            transform_op_cvec(:,:,i)=lat_vec(i,1)*op(:,:,1) + lat_vec(i,2)*op(:,:,2) + lat_vec(i,3)*op(:,:,3) 
        end forall

    end function

    function dot_prod_operator(op1,op2,lat_vec)
        complex(dpc), dimension(:,:,:),                         intent(in)      :: op1, op2
        real(dp), dimension(3,3),optional,                      intent(in)      :: lat_vec

        complex(dpc), dimension(size(op1,2),size(op2,2))                        :: dot_prod_operator


        complex(dpc),dimension(size(op1,1),size(op1,2),3)                       :: op1_cart
        complex(dpc),dimension(size(op2,1),size(op2,2),3)                       :: op2_cart
                                        
        if (size(op1,1)/=size(op2,1)) then; write(6,*) "operator mismatch dim"; stop; endif 
        

        if (present(lat_vec)) then
            write(6,*) '==<', size(op1)
            op1_cart=transform(op1,lat_vec)
            op2_cart=transform(op2,lat_vec)
            
            
            dot_prod_operator=  matmul(conjg(transpose(op1_cart(:,:,1))),op2_cart(:,:,1)) &
                               +matmul(conjg(transpose(op1_cart(:,:,2))),op2_cart(:,:,2)) &
                               +matmul(conjg(transpose(op1_cart(:,:,3))),op2_cart(:,:,3))

        else
            dot_prod_operator=  matmul(conjg(transpose(op1(:,:,1))),op2(:,:,1)) &
                               +matmul(conjg(transpose(op1(:,:,2))),op2(:,:,2)) &
                               +matmul(conjg(transpose(op1(:,:,3))),op2(:,:,3))
        endif

    end function

!    function p_block(matrices, iblock,jblock, p)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, jblock
!
!        complex(dpc), dimension(:,:,:), pointer                                 :: p_block
!
!
!        allocate(p_block(matrices%block(iblock)%band_range(2)-   &
!                   matrices%block(iblock)%band_range(1)+1, &
!                   matrices%block(jblock)%band_range(2)-   &
!                   matrices%block(jblock)%band_range(1)+1,3))
!
!        p_block=matrices%p(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                     matrices%block(jblock)%band_range(1):matrices%block(jblock)%band_range(2),:)
!        
!    
!
!    end function
!
!
!    function pk_block(matrices, iblock,jblock, q, pk)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, jblock, q
!        
!        complex(dpc), dimension(:,:,:), pointer                                 :: pk_block
!
!
!        allocate(pk_block(matrices%block(iblock)%band_range(2)-   &
!                    matrices%block(iblock)%band_range(1)+1, &
!                    matrices%dblock(2,jblock,q)-   &
!                    matrices%dblock(1,jblock,q)+1,3))
!       
!
!        pk_block=matrices%p_kdiff(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                            matrices%dblock(1,jblock,q):matrices%dblock(2,jblock,q),:,q)
!        
!    
!
!    end function
!
!
!    function sk_block(matrices, iblock, q, sk)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, q
!        
!        complex(dpc), dimension(:,:),                 pointer                   :: sk_block
!
!        allocate(sk_block(matrices%block(iblock)%band_range(2)-   &
!                    matrices%block(iblock)%band_range(1)+1, &
!                    matrices%dblock(2,iblock,q)-   &
!                    matrices%dblock(1,iblock,q)+1))
!
!        sk_block=matrices%s_kdiff(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                            matrices%dblock(1,iblock,q):matrices%dblock(2,iblock,q),q)
!        
!    
!
!    end function
!
!
!

!
!    function p_block(matrices, iblock,jblock)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, jblock
!
!        complex(dpc), dimension(matrices%block(iblock)%band_range(2)-   &
!                   matrices%block(iblock)%band_range(1)+1, &
!                   matrices%block(jblock)%band_range(2)-   &
!                   matrices%block(jblock)%band_range(1)+1,3)                    :: p_block
!
!
!
!        p_block=matrices%p(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                     matrices%block(jblock)%band_range(1):matrices%block(jblock)%band_range(2),:)
!        
!    
!
!    end function
!
!
!    function pk_block(matrices, iblock,jblock, q)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, jblock, q
!        
!        complex(dpc), dimension(matrices%block(iblock)%band_range(2)-   &
!                    matrices%block(iblock)%band_range(1)+1, &
!                    matrices%dblock(2,jblock,q)-   &
!                    matrices%dblock(1,jblock,q)+1,3)                            :: pk_block
!
!
!       
!
!        pk_block=matrices%p_kdiff(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                            matrices%dblock(1,jblock,q):matrices%dblock(2,jblock,q),:,q)
!        
!    
!
!    end function
!
!
!    function sk_block(matrices, iblock, q)
!        type(operators),                                        intent(in)      :: matrices
!        integer,                                                intent(in)      :: iblock, q
!        
!        complex(dpc), dimension(matrices%block(iblock)%band_range(2)-   &
!                    matrices%block(iblock)%band_range(1)+1, &
!                    matrices%dblock(2,iblock,q)-   &
!                    matrices%dblock(1,iblock,q)+1)                              :: sk_block
!
!
!        sk_block=matrices%s_kdiff(matrices%block(iblock)%band_range(1):matrices%block(iblock)%band_range(2), &
!                            matrices%dblock(1,iblock,q):matrices%dblock(2,iblock,q),q)
!        
!    
!
!    end function

     subroutine dump_op_kl(matrices,ispin,cutoff)
        USE mp_global, ONLY : proc_id
        type(kpt_ops),intent(in)                                                :: matrices
        integer,intent(in)                                                      :: ispin
        real(dp),intent(in)                                                     :: cutoff 

        integer                                                                 :: dumpunit, st, iblock, jblock
        character(LEN=3)                                                        :: procstr
        character(LEN=256), dimension(3)                                        :: dump_file
        integer, dimension(3)                                                   :: vshape
        integer                                                                 :: ii,jj,kk

        dump_file(1) = 'dump'
        dump_file(2) = 'dump_up'
        dump_file(3) = 'dump_dw'
        dumpunit=20+proc_id

        write(procstr,'(I3.3)')proc_id
        open(dumpunit, file=TRIM(dump_file(ispin))//procstr, form='formatted',position='append')


        write(dumpunit, '("kpt = ",3F10.5)') matrices%k
        do iblock=1,matrices%nblocks
            do jblock=1,matrices%nblocks
                if ( all(matrices%block(iblock,jblock)%band_range /= 0) &
                  .and. abs(matrices%block(iblock,jblock)%energy) < cutoff &
                  .and. abs(matrices%block(iblock,jblock)%filling) > 0d0 ) then
                    
                    write(dumpunit,'(6I2)') matrices%block(iblock,jblock)%band_range, iblock, jblock
                    write(dumpunit,'(2E15.6)') matrices%block(iblock,jblock)%filling, &
                                               matrices%block(iblock,jblock)%energy
                    vshape = shape(matrices%block(iblock,jblock)%p%v)
                    do kk=1,vshape(3)
                      do jj=1,vshape(2)
                        do ii=1,vshape(1)
                          write(dumpunit,'(2E15.6)') real(matrices%block(iblock,jblock)%p%v(ii,jj,kk)), &
                                                     aimag(matrices%block(iblock,jblock)%p%v(ii,jj,kk))
                        enddo
                      enddo
                    enddo
                endif 
            enddo
        enddo
     end subroutine


     SUBROUTINE release_matrices_data(matrices)
         type(kpt_ops)                                                        :: matrices

         integer                                                              :: iblock, jblock, is, ip

         if ( allocated(matrices%block) ) then
 
             do jblock = 1, matrices%nblocks
                 do is = 1, 3
                     deallocate(matrices%block(jblock,jblock)%s_dk(is)%s)
                     deallocate(matrices%block(jblock,jblock)%s_dk(-is)%s)
                 end do
             deallocate(matrices%block(jblock,jblock)%s_dk)
             do iblock = 1, matrices%nblocks
                 if ( associated(matrices%block(iblock,jblock)%p) )   then
                     deallocate(matrices%block(iblock,jblock)%p%v)
                     deallocate(matrices%block(iblock,jblock)%p)
                     if ( associated(matrices%block(iblock,jblock)%p_dk) ) then
                         do ip = 1, 3
                             if (allocated(matrices%block(iblock,jblock)%p_dk(ip)%v)) &
                                deallocate(matrices%block(iblock,jblock)%p_dk(ip)%v)
                             if (allocated(matrices%block(iblock,jblock)%p_dk(ip)%v)) &
                                deallocate(matrices%block(iblock,jblock)%p_dk(-ip)%v)
                         end do
                         deallocate(matrices%block(iblock,jblock)%p_dk)
                     end if
                 end if
             end do
             end do

             deallocate(matrices%block)

         end if

     END SUBROUTINE


end module
