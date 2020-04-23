module es_tools
use flinal_tools
  USE constants
  ! --------------------------------------------------------------------------
implicit none

  type lattice_data
      real(dp), dimension(3,3)                                           :: rprimd,gprimd
      real(dp)                                                           :: vol
      integer                                                            :: nsym
      integer,            allocatable, dimension(:,:,:)                  :: rsym,gsym
      real(dp),           allocatable, dimension(:,:)                    :: lsym
      complex(dpc),       allocatable, dimension(:,:,:)                  :: ssym
  end type 


  type pwset
    real(dp), dimension(3)                                             :: k
    integer,            allocatable, dimension(:,:)                    :: vec
    integer,            allocatable, dimension(:,:,:)                  :: ind
  end type

  type wavefunction
    complex(dpc),       allocatable, dimension(:,:)                    :: coeff
    type(pwset),        pointer                                        :: g    
  end type


  type block
    type(wavefunction),    allocatable, dimension(:)                   :: wf
    integer, dimension(2)                                              :: band_range
    real(dp)                                                           :: filling
    real(dp)                                                           :: energy
  end type



  type kpt_data
    type(block), dimension(:), allocatable                             :: block
    integer,                         dimension(3)                      :: k
    integer,                         dimension(-3:3)                   :: kdiff
    integer                                                            :: nband,vband,nblock, vblock
    integer                                                            :: npw,nspinor
    type(pwset),        pointer                                        :: g
  end type


  type es_data
    integer                                                            :: nkpt
    type(kpt_data),     allocatable, dimension(:)                      :: kpt
    integer,            allocatable, dimension(:,:,:)                  :: kpt_ind
    real(dp),           allocatable, dimension(:,:,:)                  :: kptns
    integer,            allocatable, dimension(:,:)                    :: rk_id
    integer,                         dimension(3)                      :: kptdim, fftdim
    real(dp),                        dimension(3)                      :: kptoffset
    type(lattice_data)                                                 :: lattice
    real(dp)                                                           :: freq_cutoff
  end type    
  

  type es_header
    integer                                                            :: nkpt
    integer                                                            :: nspin
    integer                                                            :: nspinor
    integer                                                            :: nband
    integer                                                            :: nelec
    integer,            allocatable, dimension(:,:,:)                  :: kpt_ind
    integer,            allocatable, dimension(:,:)                    :: kgrid
    integer,            allocatable, dimension(:,:)                    :: kpt_id_list
    integer,            allocatable, dimension(:,:,:)                  :: kpt_id_link
    real(dp),           allocatable, dimension(:,:,:)                  :: kptns
    real(dp),           allocatable, dimension(:,:,:)                  :: energy,filling
    integer,                         dimension(3)                      :: kptdim, fftdim
    real(dp),                        dimension(3)                      :: kptoffset
    type(kpt_data),                  dimension(-3:3)                   :: kpt
    type(lattice_data)                                                 :: lattice
    real(dp)                                                           :: freq_cutoff
    logical                                                            :: gen_time_reversal
    logical                                                            :: spinorb
    character(len=256)                                                 :: kpts_type
  end type
    







integer :: krange(3),np_grid(3),kstart,kend


interface g_
    module procedure copy_g_ref,new_g_ref
end interface

interface g
    module procedure copy_g,new_g
end interface

interface new_wf 
    module procedure new_wf_ref
end interface
    
interface kpt
    module procedure kpt_copy, kpt_sym
end interface



contains

pure function init_pauli()
    complex(dpc), dimension(2,2,4)                                           :: init_pauli
    

    init_pauli(:,1,4)=(/dcmplx(1d0,0d0)             , dcmplx(0d0,0d0)            /)
    init_pauli(:,2,4)=(/dcmplx(0d0,0d0)             , dcmplx(1d0,0d0)            /)

    init_pauli(:,1,3)=(/dcmplx(1d0,0d0)             , dcmplx(0d0,0d0)            /)
    init_pauli(:,2,3)=(/dcmplx(0d0,0d0)             , dcmplx(-1d0,0d0)           /)

    init_pauli(:,1,2)=(/dcmplx(0d0,0d0)             , dcmplx(0d0,1d0)            /)
    init_pauli(:,2,2)=(/dcmplx(0d0,-1d0)            , dcmplx(0d0,0d0)            /)

    init_pauli(:,1,1)=(/dcmplx(0d0,0d0)             , dcmplx(1d0,0d0)            /)
    init_pauli(:,2,1)=(/dcmplx(1d0,0d0)             , dcmplx(0d0,0d0)            /)

    

end function

pure function O3toSU2(o3)
    complex(dpc), dimension(2,2)                                                    :: O3toSU2

    real(dp), dimension(3,3),                                       intent(in)      :: o3
    real(dp), dimension(3,3)                                                        :: so3

    complex(dpc), dimension(2,2,4)                                                  :: pauli

    pauli=init_pauli()

    IF ( o3(1,1)*(o3(2,2)*o3(3,3)-o3(2,3)*o3(3,2)) &
       + o3(1,2)*(o3(2,3)*o3(3,1)-o3(2,1)*o3(3,3)) &
       + o3(1,3)*(o3(2,1)*o3(3,2)-o3(2,2)*o3(3,1)) < 0 ) then
        so3=-o3
    else
        so3=o3
    endif

    O3toSU2=sqrt( max( 0d0, 1d0 + so3(1,1) + so3(2,2) + so3(3,3) ) ) / 2 * pauli(:,:,4) &
            +dcmplx(0,1d0)*( sign(sqrt( max( 0d0, 1d0 + so3(1,1) - so3(2,2) - so3(3,3) ) ) / 2,so3(3,2)-so3(2,3)) * pauli(:,:,1) &  
                            +sign(sqrt( max( 0d0, 1d0 - so3(1,1) + so3(2,2) - so3(3,3) ) ) / 2,so3(1,3)-so3(3,1)) * pauli(:,:,2) &  
                            +sign(sqrt( max( 0d0, 1d0 - so3(1,1) - so3(2,2) + so3(3,3) ) ) / 2,so3(2,1)-so3(1,2)) * pauli(:,:,3) )   

end function

!!lattice functions!!!!!!
    function lattice(rprimd, rsym, lsym)
        type(lattice_data)                                                      :: lattice

        real(dp), dimension(3,3),                               intent(in)      :: rprimd
        integer, dimension(:,:,:),                              intent(in)      :: rsym
        real(dp), dimension(:,:),                               intent(in)      :: lsym

        integer                                                                 :: isym

        !lattice shape/size
        lattice%rprimd=rprimd
        lattice%gprimd=transpose(rinverse(rprimd))
        lattice%vol=abs(det(rprimd))

        !symmetries
        lattice%nsym=size(rsym,3)
        ALLOCATE( lattice%rsym(3,3,lattice%nsym), lattice%gsym(3,3,lattice%nsym), &
                  lattice%lsym(3,lattice%nsym), lattice%ssym(2,2,lattice%nsym) )
        lattice%rsym = rsym
        lattice%lsym = lsym
        do isym=1,lattice%nsym
            lattice%gsym(:,:,isym)=transpose(real(rsym(:,:,isym), kind=dp))
    !       lattice%ssym(:,:,isym)=O3toSU2( matmul(transpose(lattice%rprimd),matmul(rsym(:,:,isym),lattice%gprimd)))
            lattice%ssym(:,:,isym)=O3toSU2( matmul(         (lattice%rprimd),matmul(rsym(:,:,isym),transpose(lattice%gprimd))))
        end do

    end function

    subroutine lattice_bcast(lattice, root, comm)
        USE mp_global, ONLY : proc_id, mp_barrier, mp_bcast
        type(lattice_data),                                     intent(inout)   :: lattice
        integer,                                                intent(in)      :: root, comm
    
        CALL mp_barrier(comm)

        CALL mp_bcast(lattice%rprimd, root, comm)
        CALL mp_bcast(lattice%gprimd, root, comm)
        CALL mp_bcast(lattice%vol, root, comm)

        CALL mp_bcast(lattice%nsym, root, comm)

        if (proc_id/= root) then
            allocate(lattice%gsym(3,3,lattice%nsym))
            allocate(lattice%rsym(3,3,lattice%nsym))
            allocate(lattice%lsym(3,lattice%nsym))
            allocate(lattice%ssym(2,2,lattice%nsym))
        endif

        CALL mp_bcast(lattice%gsym, root, comm)
        CALL mp_bcast(lattice%rsym, root, comm)
        CALL mp_bcast(lattice%lsym, root, comm)
        CALL mp_bcast(lattice%ssym, root, comm)
    end subroutine


!!kpt functions!!!!!!
  SUBROUTINE kpt_init_tr( kcoord, kpt, gvec, bands, energy, filling, tr )
    !
    IMPLICIT NONE
    !
        type(kpt_data),                                intent(inout)       :: kpt

        integer, dimension(:,:)            ,           intent(in)          :: gvec
        complex(dpc), dimension(:,:,:),                intent(in)          :: bands
        real(dp), dimension(:),                        intent(in)          :: energy,filling
        real(dp), dimension(3),                        intent(in)          :: kcoord
        logical,                                       intent(in)          :: tr


!       complex(dpc), dimension(size(bands,1)*2,size(bands,2),size(bands,3))           :: bands_tmp
        complex(dpc), dimension(:,:,:),allocatable           :: bands_tmp
        integer, dimension(3,size(gvec,2)*2)                                :: gvec_tmp
        
        integer                                                             :: ipw, nbnd_tmp, ibnd,ispin,npw_tmp


        if (tr) then

            allocate( bands_tmp(size(bands,1)*2,size(bands,2),size(bands,3)) )
            npw_tmp = size(bands,1)
            bands_tmp(:npw_tmp,:,:)=bands(:npw_tmp,:,:)
            bands_tmp(npw_tmp+1:,:,:)=conjg(bands(:npw_tmp,:,:))

            gvec_tmp(:,:size(gvec,2))=gvec
            forall(ipw=1:size(gvec,2)) gvec_tmp(:,size(gvec,2)+ipw)=-gvec(:,ipw)-nint(2*kcoord)
           
            if (all(nint(2*kcoord)==0)) then
!                if (present(cutoff)) then
                    call kpt_init(kcoord,kpt,gvec_tmp(:,2:),bands_tmp(2:,:,:),energy,filling)
!                else
!                    call kpt_init(kcoord,kpt,gvec_tmp(:,2:),bands_tmp(2:,:,:),energy,filling,block_tol,sg)
!                endif
            else
!                if (present(cutoff)) then
                    call kpt_init(kcoord,kpt,gvec_tmp,bands_tmp,energy,filling)
!                else
!                    call kpt_init(kcoord,kpt,gvec_tmp,bands_tmp,energy,filling,block_tol,sg)
!                endif
            endif 

        else
            
!            if (present(cutoff)) then
                call kpt_init(kcoord,kpt,gvec,bands,energy,filling)
!            else
!                call kpt_init(kcoord,kpt,gvec,bands,energy,filling,block_tol,sg)
!            endif

        endif
        if (allocated(bands_tmp)) deallocate( bands_tmp )


    end subroutine

    !###modified in 4/22/12
  SUBROUTINE kpt_init( kcoord, kpt, gvec, bands, energy, filling )
    USE control, ONLY : block_tol, cutoff, sg_calc
        type(kpt_data),                                intent(inout)       :: kpt

        integer, dimension(:,:),                       intent(in)          :: gvec
        complex(dpc), dimension(:,:,:),                intent(in)          :: bands
        real(dp), dimension(:),                        intent(in)          :: energy,filling
        real(dp), dimension(3),                        intent(in)          :: kcoord

        type(block), dimension(:),      allocatable                        :: block_temp
        integer                                                            :: ipw, ibnd, vband, lvband, ucband, iblock, vblock

        !****
        integer                                                            :: ndbnd, jbnd
        !****
        if (allocated(kpt%block)) deallocate(kpt%block)
        kpt%npw=size(bands,1)
        kpt%nspinor=size(bands,2)
        kpt%nband=size(bands,3)

        if (associated(kpt%g)) then
             deallocate(kpt%g%vec)
             deallocate(kpt%g%ind)
             deallocate(kpt%g)
        endif

        kpt%g=>g_(kcoord,gvec) 

        kpt%vband=count(filling > eq_tol)
        if (.not. sg_calc) then
!            if (present(cutoff)) then
                if (energy(kpt%nband)-energy(kpt%vband)<cutoff) then
        !           write(6,*) "not enough bands for cutoff", kpt%nband, kpt%vband, energy(kpt%nband), energy(kpt%vband), cutoff
        !           call flush(6)
        !           stop
                endif
                lvband=count(energy(kpt%vband+1)-energy>cutoff)+1
!            else
!                lvband=1
!            endif
        else
!            if (present(cutoff)) then
                if (energy(kpt%nband)-energy(kpt%vband)<cutoff) then
        !           write(6,*) "not enough bands for cutoff", kpt%nband, kpt%vband, energy(kpt%nband), energy(kpt%vband), cutoff
        !           call flush(6)
        !           stop
                endif
!            end if
            lvband=1
        end if

        allocate(block_temp(kpt%nband))
        !****
        ndbnd=0
        iblock=1
        forall(ibnd=1:size(energy))
            block_temp(ibnd)%band_range=0
            block_temp(ibnd)%filling=0d0
            block_temp(ibnd)%energy=0d0
        end forall

        do ibnd=lvband,size(energy)

            if ( (ibnd<size(energy))  .and.(abs(energy(ibnd)-energy(ibnd+1))    <block_tol) ) then
                ndbnd=ndbnd+1
            else
                block_temp(iblock)%band_range=(/ ibnd-ndbnd,ibnd /)
                block_temp(iblock)%filling=sum(filling(ibnd-ndbnd:ibnd))/(1d0*ndbnd+1d0)
                block_temp(iblock)%energy=sum(energy(ibnd-ndbnd:ibnd))/(1d0*ndbnd+1d0)

                if (kpt%vband <= ibnd .and. kpt%vband >= ibnd-ndbnd) kpt%vblock=iblock

                allocate(block_temp(iblock)%wf(ibnd-ndbnd:ibnd))
                do jbnd=ibnd-ndbnd,ibnd
                    allocate(block_temp(iblock)%wf(jbnd)%coeff(size(bands,1),size(bands,2)) )
                    block_temp(iblock)%wf(jbnd)%coeff=bands(:,:,jbnd)!new_wf(bands(:,:,jbnd),g)
                enddo

                iblock=iblock+1
                ndbnd=0
            endif
        enddo
        !****
      ! block_temp= blocks(bands,        energy,filling,lvband,block_tol)
        
        kpt%nblock=count(block_temp%band_range(1)>0)
        allocate(kpt%block(kpt%nblock))
        kpt%block=block_temp(1:kpt%nblock)

        do iblock=1,size(kpt%block)
            do ibnd=kpt%block(iblock)%band_range(1),kpt%block(iblock)%band_range(2)
                kpt%block(iblock)%wf(ibnd)%g=>kpt%g
            end do
        end do





        if (allocated(block_temp)) then
            do iblock=1,size(block_temp)
                if (allocated(block_temp(iblock)%wf)) then
                    do ibnd=block_temp(iblock)%band_range(1),block_temp(iblock)%band_range(2)
                        if (allocated(block_temp(iblock)%wf(ibnd)%coeff))  deallocate(block_temp(iblock)%wf(ibnd)%coeff)
                        if (associated(block_temp(iblock)%wf(ibnd)%g))     deallocate(block_temp(iblock)%wf(ibnd)%g)
                    enddo
                endif
            enddo                
            deallocate(block_temp)
        endif


    end subroutine

    subroutine release_kpt(kpt)
        type(kpt_data),                                     intent(inout)       :: kpt
        

        integer                                                                 :: iblock, ibnd

        if (allocated(kpt%block)) then
            do iblock=1,kpt%nblock
                if (allocated(kpt%block(iblock)%wf)) then
                    do ibnd=kpt%block(iblock)%band_range(1),kpt%block(iblock)%band_range(2)
                        if (allocated(kpt%block(iblock)%wf(ibnd)%coeff))   deallocate(kpt%block(iblock)%wf(ibnd)%coeff)
            !           if (associated(kpt%block(iblock)%wf(ibnd)%g))     deallocate(kpt%block(iblock)%wf(ibnd)%g)
                    enddo
                endif
            enddo                
            deallocate(kpt%block)
        endif
        if (associated(kpt%g)) deallocate(kpt%g)
    end subroutine

    !###added in 4/22/12
    SUBROUTINE release_kpt_kl(es_kl)
        type(es_header),                                    intent(inout)       :: es_kl
        

        integer                                                                 :: iblock, ibnd, ikpt

        do ikpt=-3,3
            if (allocated(es_kl%kpt(ikpt)%block)) then
                do iblock=1,es_kl%kpt(ikpt)%nblock
                    if (allocated(es_kl%kpt(ikpt)%block(iblock)%wf)) then
                        do ibnd=es_kl%kpt(ikpt)%block(iblock)%band_range(1),es_kl%kpt(ikpt)%block(iblock)%band_range(2)
                            IF (allocated(es_kl%kpt(ikpt)%block(iblock)%wf(ibnd)%coeff)) &
                              deallocate(es_kl%kpt(ikpt)%block(iblock)%wf(ibnd)%coeff)
                        enddo
                    endif
                enddo                
                deallocate(es_kl%kpt(ikpt)%block)
            endif
            if (associated(es_kl%kpt(ikpt)%g)) deallocate(es_kl%kpt(ikpt)%g)
        end do

    end subroutine


    subroutine kpt_bcast(kpt, root, comm)
        USE mp_global, ONLY : proc_id, mp_bcast
        type(kpt_data),                                         intent(inout)   :: kpt
        integer,                                                intent(in)      :: root, comm

        integer, dimension(2,3)                                                 :: ginddim
        integer                                                                 :: i, ibnd, iblock
        logical                                                                 :: band_is_stored

        CALL mp_bcast(kpt%k, root, comm)
!        call MPI_BCAST(kpt%knext,3,MPI_INTEGER, root, comm, ierr)
!        call MPI_BCAST(kpt%kprev,3,MPI_INTEGER, root, comm, ierr)
        
        CALL mp_bcast(kpt%npw, root, comm)
        CALL mp_bcast(kpt%nspinor, root, comm)
        CALL mp_bcast(kpt%nband, root, comm)
        CALL mp_bcast(kpt%vband, root, comm)
        CALL mp_bcast(kpt%vblock, root, comm)
        CALL mp_bcast(kpt%nblock, root, comm)

        if (proc_id == root) then
            forall (i=1:3) ginddim(:,i)=(/ lbound(kpt%g%ind,i),ubound(kpt%g%ind,i) /)
        endif
        CALL mp_bcast(ginddim, root, comm)

        if (proc_id /= root) then
            if (associated(kpt%g)) deallocate(kpt%g)

            ALLOCATE( kpt%g )
            ALLOCATE( kpt%g%vec(3,kpt%npw) )
            ALLOCATE( kpt%g%ind(ginddim(1,1):ginddim(2,1), &
                                ginddim(1,2):ginddim(2,2), &
                                ginddim(1,3):ginddim(2,3)) )
        endif
        CALL mp_bcast(kpt%g%vec, root, comm)
        CALL mp_bcast(kpt%g%ind, root, comm)
        CALL mp_bcast(kpt%g%k, root, comm)


        if (proc_id /= root) then
            if (allocated(kpt%block)) deallocate(kpt%block)

            allocate(kpt%block(kpt%nblock))
        endif

        do iblock=1,kpt%nblock
            CALL mp_bcast(kpt%block(iblock)%band_range, root, comm)
            CALL mp_bcast(kpt%block(iblock)%filling, root, comm)
            CALL mp_bcast(kpt%block(iblock)%energy, root, comm)

            if (proc_id /= root) allocate(kpt%block(iblock)%wf(kpt%block(iblock)%band_range(1):kpt%block(iblock)%band_range(2)))
            do ibnd=kpt%block(iblock)%band_range(1),kpt%block(iblock)%band_range(2)
                if (proc_id /= root) allocate(kpt%block(iblock)%wf(ibnd)%coeff(kpt%npw, kpt%nspinor))
                CALL mp_bcast(kpt%block(iblock)%wf(ibnd)%coeff, root, comm)
                kpt%block(iblock)%wf(ibnd)%g=>kpt%g
            enddo
            
        enddo
            
    end subroutine

    function kpt_copy(kpt)
        type(kpt_data)                                                          :: kpt_copy

        type(kpt_data),                                  intent(in)             :: kpt  

        integer                                                                 :: ibnd, iblock 

        kpt_copy%npw=kpt%npw
        kpt_copy%nspinor=kpt%nspinor
        kpt_copy%nband=kpt%nband
        kpt_copy%vband=kpt%vband
        kpt_copy%vblock=kpt%vblock
        kpt_copy%nblock=kpt%nblock

        kpt_copy%g=>g_(kpt%g)
!!!TODO: maybe flesh out more copy constructors

        allocate(kpt_copy%block(kpt_copy%nblock))
        do iblock=1,kpt_copy%nblock
            kpt_copy%block(iblock)=block_copy(kpt%block(iblock),kpt_copy%g)
        enddo

    end function

    function kpt_sym(kpt,gsym,lsym,ssym,tr)
        type(kpt_data)                                                          :: kpt_sym

        type(kpt_data),                                 intent(inout)           :: kpt  
        integer, dimension(3,3),                        intent(in)              :: gsym
        real(dp), dimension(3),                         intent(in)              :: lsym
        complex(dpc), dimension(2,2),                   intent(in)              :: ssym
        logical,                                        intent(in)              :: tr

        integer, dimension(3,3)                                                 :: gsym_eff
        complex(dpc), dimension(:,:), allocatable                               :: lsym_factor
        integer                                                                 :: ipw, ibnd,iblock
        integer, dimension(3)                                                   :: kshift

        kpt_sym=kpt_copy(kpt)


        if (tr) gsym_eff=-gsym
        if (.not.tr) gsym_eff=gsym
       
        kshift=nint(nmod(matmul(gsym_eff,kpt_sym%g%k),1d0)+eq_tol  -matmul(gsym_eff,kpt_sym%g%k))
        
        
        kpt_sym%g%k=matmul(gsym_eff,kpt_sym%g%k)+real(kshift,kind=dp)
        forall (ipw=1:kpt_sym%npw) kpt_sym%g%vec(:,ipw)=matmul(gsym_eff,kpt_sym%g%vec(:,ipw))-kshift
         

        kpt_sym%g=g(kpt_sym%g%k,kpt_sym%g%vec)


        !!!!rotate spinor
        if (kpt_sym%nspinor==2) then
            forall (iblock=1:kpt_sym%nblock) 
                forall(ibnd=kpt_sym%block(iblock)%band_range(1):kpt_sym%block(iblock)%band_range(2)) 
                    kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:)=matmul(kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:),transpose(ssym))
                endforall
            endforall
        endif
        

        allocate(lsym_factor(kpt_sym%npw,kpt_sym%nspinor))
        forall (ipw=1:kpt_sym%npw)  lsym_factor(ipw,:)=exp(cmplx(0d0,1d0,kind=dpc)*2*PI*dot_product(lsym,kpt_sym%g%vec(:,ipw)))

        if (tr) then
            forall (iblock=1:kpt_sym%nblock) 
                forall(ibnd=kpt_sym%block(iblock)%band_range(1):kpt_sym%block(iblock)%band_range(2)) 
                    kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:)=lsym_factor*conjg(kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:))
                endforall
            endforall
        else
            forall (iblock=1:kpt_sym%nblock)
                forall (ibnd=kpt_sym%block(iblock)%band_range(1):kpt_sym%block(iblock)%band_range(2)) 
                    kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:)=lsym_factor*kpt_sym%block(iblock)%wf(ibnd)%coeff(:,:)
                endforall
            end forall
        endif
        
        call release_kpt(kpt)
        if (allocated(lsym_factor)) deallocate(lsym_factor)

    end function


    !biases the modulo(needs rewrite)
    elemental real(dp) function nmod(x,div)
        real(dp),                                           intent(in)      :: x,div

        nmod=x-nint((x-eq_tol)/div)*div+eq_tol

    end function




    logical function is_tr(kcoord)
        real(dp), dimension(3),                             intent(in)      :: kcoord

!        is_tr=all(abs(modulo(kcoord,1d0)-modulo(-kcoord,1d0)) <= eq_tol)
        is_tr=all(abs(modulo(abs(2*kcoord)+.5d0,1d0)-.5) <= eq_tol)

    end function

!!block functions!!!!!!
    function block_copy(inblock,g)
        type (block)                                                        :: block_copy
        type (block),                                       intent(in)      :: inblock

        type (pwset), pointer,                              intent(in)      :: g


        integer                                                             :: ibnd        

        block_copy%band_range=inblock%band_range
        block_copy%energy=inblock%energy
        block_copy%filling=inblock%filling

        
        allocate(block_copy%wf(inblock%band_range(1):inblock%band_range(2)))
        do ibnd=inblock%band_range(1),inblock%band_range(2)
            block_copy%wf(ibnd)=new_wf(inblock%wf(ibnd)%coeff, g)
        enddo
            


    end function


    function blocks(bands,    energy,filling,lvband,tol)
        complex(dpc), dimension(:,:,:),                     intent(in)      :: bands
!       type (pwset), pointer,                              intent(in)      :: g
        real(dp),dimension(:),                              intent(in)      :: energy,filling
        integer,                                            intent(in)      :: lvband
        real(dp),                                           intent(in)      :: tol

        type(block), dimension(size(energy))                                :: blocks
        

        integer                                                             :: ibnd,jbnd, ikpt_n,iblock,ndbnd


        ndbnd=0
        iblock=1
        forall(ibnd=1:size(energy)) 
            blocks(ibnd)%band_range=0
            blocks(ibnd)%filling=0d0
            blocks(ibnd)%energy=0d0
        end forall

        do ibnd=lvband,size(energy)
            
            if ( (ibnd<size(energy))  .and.(abs(energy(ibnd)-energy(ibnd+1))    <tol) ) then
                ndbnd=ndbnd+1
            else
                blocks(iblock)%band_range=(/ ibnd-ndbnd,ibnd /)
                blocks(iblock)%filling=sum(filling(ibnd-ndbnd:ibnd))/(1d0*ndbnd+1d0)
                blocks(iblock)%energy=sum(energy(ibnd-ndbnd:ibnd))/(1d0*ndbnd+1d0)

                allocate(blocks(iblock)%wf(ibnd-ndbnd:ibnd))
                do jbnd=ibnd-ndbnd,ibnd
                    allocate(blocks(iblock)%wf(jbnd)%coeff(size(bands,1),size(bands,2)) )
    !       !       blocks(iblock)%wf(jbnd)%coeff=bands(:,:,jbnd)!new_wf(bands(:,:,jbnd),g)
    !    !  !       blocks(iblock)%wf(jbnd)%g=>g
                enddo

                iblock=iblock+1
                ndbnd=0
            endif
        enddo
    end function


!!es_construction!!!!

  ! --------------------------------------------------------------------------
  SUBROUTINE es_header_bcast(es, root, comm)
  ! --------------------------------------------------------------------------
    USE mp_global, ONLY : proc_id, mp_bcast
    !
    IMPLICIT NONE
    !
    TYPE(es_header), POINTER, INTENT(inout) :: es
    INTEGER, INTENT(in) :: root, comm
    !
    IF (proc_id/=root) ALLOCATE(es)
    !
    CALL lattice_bcast(es%lattice, root, comm) 
    !
    CALL mp_bcast(es%nkpt, root, comm)
    CALL mp_bcast(es%nspin, root, comm)
    CALL mp_bcast(es%nspinor, root, comm)
    CALL mp_bcast(es%nband, root, comm)
    CALL mp_bcast(es%nelec, root, comm)
    CALL mp_bcast(es%kptdim, root, comm)
    CALL mp_bcast(es%kptoffset, root, comm)
    CALL mp_bcast(es%fftdim, root, comm)
    CALL mp_bcast(es%freq_cutoff, root, comm)
    CALL mp_bcast(es%gen_time_reversal, root, comm)
    CALL mp_bcast(es%spinorb, root, comm)
    !
    CALL mp_bcast(es%kpts_type, root, comm)
    IF (proc_id/=root) ALLOCATE ( es%kptns(3, es%lattice%nsym, es%nkpt) )
    CALL mp_bcast(es%kptns, root, comm)
    !
  END SUBROUTINE es_header_bcast
  !
  ! --------------------------------------------------------------------------
  SUBROUTINE energy_filling_bcast( es, root, comm )
  ! --------------------------------------------------------------------------
    USE mp_global, ONLY : proc_id, mp_bcast
    !
    IMPLICIT NONE
    !
    TYPE(es_header), INTENT(inout) :: es
    INTEGER, INTENT(in) :: root, comm
    !
    ! Allocate the arrays
    IF (proc_id/=root) THEN
      ALLOCATE( es%energy(es%nband, es%nspin, es%nkpt), &
                es%filling(es%nband, es%nspin, es%nkpt) )
    END IF
    !
    CALL mp_bcast(es%energy, root, comm)
    CALL mp_bcast(es%filling, root, comm)
    !
  END SUBROUTINE energy_filling_bcast



            
!   subroutine es_symkpts(es,kpt,tr)
!       type(es_data),pointer,                              intent(inout)       :: es
!       type(kpt_data),                                     intent(in)          :: kpt
!       logical,                                            intent(in)          :: tr

!       real(dp), dimension(3)                                                  :: kcoord
!       integer                                                                 :: k_id, k_tr
!       
!       integer                                                                 :: isym
!       
!       
!       !!TODO:kpt array should really be a linked mesh
!       do isym=1, es%lattice%nsym
!           kcoord=matmul(es%lattice%gsym(:,:,isym),kpt%g%k)
!           k_id=get_kpt_id(kcoord,es%kptoffset,es%kptdim)
!           k_tr=get_kpt_id(-kcoord,es%kptoffset,es%kptdim)

!           if (proc_id==0) write(6,'(A, I, A, 3F8.4)') "sym", isym, "  kcoord", kcoord
!           call flush(6)
!      
!           if (         any(get_kprocs(es,k_id)==proc_id) .and..not.allocated(es%kpt(k_id)%block ))           es%kpt(k_id)= kpt_sym(kpt,es%lattice%gsym(:,:,isym),es%lattice%lsym(:,isym),es%lattice%ssym(:,:,isym),.false.)
!           if (tr .and. any(get_kprocs(es,k_tr)==proc_id) .and..not.allocated(es%kpt(k_tr)%block ))           es%kpt(k_tr)= kpt_sym(kpt,es%lattice%gsym(:,:,isym),es%lattice%lsym(:,isym),es%lattice%ssym(:,:,isym),.true.)

!       enddo

!   end subroutine


        
    !###added in 4/22/12
    SUBROUTINE get_kpt_id_list(es)
      type(es_header),                               intent(inout)       :: es

      real(dp), dimension(3)                                             :: kcoord
      integer                                                            :: ikpt, isym


      if (es%gen_time_reversal) then 
          
          allocate( es%kpt_id_list(es%nkpt,2*es%lattice%nsym) ) 
          do ikpt=1,es%nkpt
              do isym=1,es%lattice%nsym
                  kcoord=matmul( es%lattice%gsym(:,:,isym),es%kptns(:,1,ikpt) )
                  es%kptns(:,isym,ikpt)=kcoord
                  es%kpt_id_list(ikpt,isym)=get_kpt_id(kcoord,es%kptoffset,es%kptdim)
                  es%kpt_id_list(ikpt,isym+es%lattice%nsym)=get_kpt_id(-kcoord,es%kptoffset,es%kptdim)
              end do
          end do
      else
          allocate( es%kpt_id_list(es%nkpt,  es%lattice%nsym) ) 
          do ikpt=1,es%nkpt
              do isym=1,es%lattice%nsym
                  kcoord=matmul( es%lattice%gsym(:,:,isym),es%kptns(:,1,ikpt) )
                  es%kptns(:,isym,ikpt)=kcoord
                  if (es%kpts_type=='nonmp' .or. es%kpts_type=='path') then
                      !ignore symmetries of non mp grid
                      es%kpt_id_list(ikpt,isym)= ikpt
                  else
                      es%kpt_id_list(ikpt,isym)=get_kpt_id(kcoord,es%kptoffset,es%kptdim)
                  end if

              end do
          end do
      end if

    end SUBROUTINE



    !###added in 4/22/12
    FUNCTION check_kpt_id_list( kpt_id_list, id)
      integer,dimension(2)                                               :: check_kpt_id_list
      integer,dimension(:,:),                       intent(in)           :: kpt_id_list
      integer,                                      intent(in)           :: id

      integer,dimension(size(kpt_id_list,1),size(kpt_id_list,2))         :: kpt_id_array_tmp

      kpt_id_array_tmp=id
      check_kpt_id_list=minloc(abs(kpt_id_list-kpt_id_array_tmp))

    end FUNCTION




    !###added in 4/22/12
    SUBROUTINE get_kpt_id_link(es,nproc)
      type(es_header),                               intent(inout)       :: es
      integer,                                       intent(in)          :: nproc

      integer                                                            :: kpt_perset, kpt_nset, ikpt
      real(dp)                                                           :: nproc_t


      nproc_t=real(nproc,kind=dp)
      kpt_perset=nproc
      if (es%kpts_type=='path') then
          kpt_nset=ceiling( es%nkpt/nproc_t )
      else
          kpt_nset=ceiling( product(es%kptdim)/nproc_t )
      endif

      allocate( es%kpt_id_link(-3:3,kpt_perset,kpt_nset) )
      es%kpt_id_link=-1

      if (es%kpts_type=='path') then
          do ikpt=1,es%nkpt
              es%kpt_id_link(0,mod(ikpt-1,kpt_perset)+1,int((ikpt-1)/kpt_perset)+1)=ikpt
          enddo

      else

          do ikpt=1,product(es%kptdim)
              es%kpt_id_link(1:3,mod(ikpt-1,kpt_perset)+1,int((ikpt-1)/kpt_perset)+1)=                     &
                                  (/mod(es%kgrid(1,ikpt),es%kptdim(1))+1                                 + &
                                    es%kptdim(1)*(es%kgrid(2,ikpt)-1)                                    + &
                                    es%kptdim(1)*es%kptdim(2)*(es%kgrid(3,ikpt)-1),                        &
    
                                    es%kgrid(1,ikpt)                                                     + &
                                    mod(es%kgrid(2,ikpt),es%kptdim(2))*es%kptdim(1)                      + &
                                    es%kptdim(1)*es%kptdim(2)*(es%kgrid(3,ikpt)-1),                        &
    
                                    es%kgrid(1,ikpt)                                                     + &
                                    es%kptdim(1)*(es%kgrid(2,ikpt)-1)                                    + &
                                    mod(es%kgrid(3,ikpt),es%kptdim(3))*es%kptdim(1)*es%kptdim(2)/)
    
              es%kpt_id_link(-3:-1,mod(ikpt-1,kpt_perset)+1,int((ikpt-1)/kpt_perset)+1)=                   &
                                  (/es%kgrid(1,ikpt)                                                     + &
                                    es%kptdim(1)*(es%kgrid(2,ikpt)-1)                                    + &
                                    modulo(es%kgrid(3,ikpt)-2,es%kptdim(3))*es%kptdim(1)*es%kptdim(2),     &
    
                                    es%kgrid(1,ikpt)                                                     + &
                                    modulo(es%kgrid(2,ikpt)-2,es%kptdim(2))*es%kptdim(1)                 + &
                                    es%kptdim(1)*es%kptdim(2)*(es%kgrid(3,ikpt)-1),                        &
    
                                    modulo(es%kgrid(1,ikpt)-2,es%kptdim(1))+1                            + &
                                    es%kptdim(1)*(es%kgrid(2,ikpt)-1)                                    + &
                                    es%kptdim(1)*es%kptdim(2)*(es%kgrid(3,ikpt)-1)/)
    
              es%kpt_id_link(0,mod(ikpt-1,kpt_perset)+1,int((ikpt-1)/kpt_perset)+1)=ikpt
    
              !if non mp grid, skip points on the boundary
              if (es%kpts_type=='nonmp') then
                  if (     (es%kgrid(1,ikpt) .eq. 1) &
                      .or. (es%kgrid(1,ikpt) .eq. es%kptdim(1)) &
                      .or. (es%kgrid(2,ikpt) .eq. 1) &
                      .or. (es%kgrid(2,ikpt) .eq. es%kptdim(2)) &
                      .or. (es%kgrid(3,ikpt) .eq. 1) &
                      .or. (es%kgrid(3,ikpt) .eq. es%kptdim(3)) ) then
    
                      es%kpt_id_link(-3:3,mod(ikpt-1,kpt_perset)+1,int((ikpt-1)/kpt_perset)+1) = -1   
                  end if
              end if
          end do
      endif



       
    end SUBROUTINE



    

    !###added in 4/12/22
    SUBROUTINE get_kpt_grid_kl(es_kl)
      type(es_header),                               intent(inout)       :: es_kl

      integer                                                            :: ikpt, i


      allocate( es_kl%kpt_ind(es_kl%kptdim(1),es_kl%kptdim(2),es_kl%kptdim(3)), es_kl%kgrid(3,product(es_kl%kptdim)) )
      do ikpt=1,product(es_kl%kptdim)
          es_kl%kgrid(:,ikpt)=(/mod(ikpt-1,es_kl%kptdim(1))+1,                                             &
                                mod((ikpt-1)/es_kl%kptdim(1),es_kl%kptdim(2))+1,                           &
                                mod((ikpt-1)/(es_kl%kptdim(1)*es_kl%kptdim(2)),es_kl%kptdim(3))+1/)
          es_kl%kpt_ind(es_kl%kgrid(1,ikpt),es_kl%kgrid(2,ikpt),es_kl%kgrid(3,ikpt))=ikpt
      end do

    end SUBROUTINE


    subroutine get_kpt_grid(es)
      type(es_data),                                 intent(inout)       :: es
   
      integer                                                            :: ikpt,i
   
   
   
      allocate(es%kpt(PRODUCT(es%kptdim)),es%kpt_ind(es%kptdim(1),es%kptdim(2),es%kptdim(3)))
      do ikpt=1,size(es%kpt)
        es%kpt(ikpt)%k=(/mod(ikpt-1,es%kptdim(1))+1,                                             &
                         mod((ikpt-1)/es%kptdim(1),es%kptdim(2))+1,                              &
                         mod((ikpt-1)/(es%kptdim(1)*es%kptdim(2)),es%kptdim(3))+1/)
   
        es%kpt_ind(es%kpt(ikpt)%k(1),es%kpt(ikpt)%k(2),es%kpt(ikpt)%k(3))=ikpt
      enddo
   
    end subroutine

    subroutine link_kpts(es)
    type(es_data),                                 intent(inout)       :: es

    integer                                                            :: ikpt

!   es%nkpt=size(es%kpt)

    do ikpt=1,size(es%kpt)
      es%kpt(ikpt)%k=(/mod(ikpt-1,es%kptdim(1))+1,                                             &
                       mod((ikpt-1)/es%kptdim(1),es%kptdim(2))+1,                              &
                       mod((ikpt-1)/(es%kptdim(1)*es%kptdim(2)),es%kptdim(3))+1/)


      es%kpt(ikpt)%kdiff(1:3)=(/mod(es%kpt(ikpt)%k(1),es%kptdim(1))+1                             + &
                                es%kptdim(1)*(es%kpt(ikpt)%k(2)-1)                                + &
                                es%kptdim(1)*es%kptdim(2)*(es%kpt(ikpt)%k(3)-1),                    &

                                es%kpt(ikpt)%k(1)                                                 + &
                                mod(es%kpt(ikpt)%k(2),es%kptdim(2))*es%kptdim(1)                  + &
                                es%kptdim(1)*es%kptdim(2)*(es%kpt(ikpt)%k(3)-1),                    &

                                es%kpt(ikpt)%k(1)                                                 + &
                                es%kptdim(1)*(es%kpt(ikpt)%k(2)-1)                                + &
                                mod(es%kpt(ikpt)%k(3),es%kptdim(3))*es%kptdim(1)*es%kptdim(2)/)

      es%kpt(ikpt)%kdiff(-3:-1)=(/es%kpt(ikpt)%k(1)                                                     + &
                                  es%kptdim(1)*(es%kpt(ikpt)%k(2)-1)                                    + &
                                  modulo(es%kpt(ikpt)%k(3)-2,es%kptdim(3))*es%kptdim(1)*es%kptdim(2),     &

                                  es%kpt(ikpt)%k(1)                                                 + &
                                  modulo(es%kpt(ikpt)%k(2)-2,es%kptdim(2))*es%kptdim(1)             + &
                                  es%kptdim(1)*es%kptdim(2)*(es%kpt(ikpt)%k(3)-1),                    &

                                  modulo(es%kpt(ikpt)%k(1)-2,es%kptdim(1))+1                        + &
                                  es%kptdim(1)*(es%kpt(ikpt)%k(2)-1)                                + &
                                  es%kptdim(1)*es%kptdim(2)*(es%kpt(ikpt)%k(3)-1)/)
      es%kpt(ikpt)%kdiff(0)=ikpt
    enddo


  end subroutine




    integer function get_kpt_id(kcoord,offset,kptdim)
        real(dp), dimension(3),                          intent(in)          :: kcoord,offset
        integer, dimension(3),                           intent(in)          :: kptdim


        get_kpt_id=1 +                     MODULO(NINT(kptdim(1)*(kcoord(1)-offset(1))),kptdim(1)) & 
                     +           kptdim(1)*MODULO(NINT(kptdim(2)*(kcoord(2)-offset(2))),kptdim(2)) &
                     + kptdim(2)*kptdim(1)*MODULO(NINT(kptdim(3)*(kcoord(3)-offset(3))),kptdim(3))

    end function

    function get_kprocs(es,ikpt)
        type(es_data),                                          intent(in)   :: es
        integer, dimension(7)                                                :: get_kprocs

        integer,                                                intent(in)   :: ikpt


        get_kprocs(1)=      (es%kpt(ikpt)%k(3)-1             )/krange(3)*np_grid(1)*np_grid(2)+  &
                            (es%kpt(ikpt)%k(2)-1             )/krange(2)*np_grid(1)+(es%kpt(ikpt)%k(1)-1)/krange(1)
        get_kprocs(2)=modulo(es%kpt(ikpt)%k(3)-2,es%kptdim(3))/krange(3)*np_grid(1)*np_grid(2)+  &
                            (es%kpt(ikpt)%k(2)-1)/krange(2)*np_grid(1)+(es%kpt(ikpt)%k(1)-1)/krange(1)
        get_kprocs(3)=modulo(es%kpt(ikpt)%k(3)  ,es%kptdim(3))/krange(3)*np_grid(1)*np_grid(2)+  &
                            (es%kpt(ikpt)%k(2)-1)/krange(2)*np_grid(1)+(es%kpt(ikpt)%k(1)-1)/krange(1)
        get_kprocs(4)=      (es%kpt(ikpt)%k(3)-1             )/krange(3)*np_grid(1)*np_grid(2)+  &
                      modulo(es%kpt(ikpt)%k(2)-2,es%kptdim(2))/krange(2)*np_grid(1)+(es%kpt(ikpt)%k(1)-1)/krange(1)
        get_kprocs(5)=      (es%kpt(ikpt)%k(3)-1             )/krange(3)*np_grid(1)*np_grid(2)+  &
                      modulo(es%kpt(ikpt)%k(2)  ,es%kptdim(2))/krange(2)*np_grid(1)+(es%kpt(ikpt)%k(1)-1)/krange(1)

        if (get_kprocs(2)==get_kprocs(1)) get_kprocs(2)=-1
        if (get_kprocs(3)==get_kprocs(1)) get_kprocs(3)=-1
        if (get_kprocs(4)==get_kprocs(1)) get_kprocs(4)=-1
        if (get_kprocs(5)==get_kprocs(1)) get_kprocs(5)=-1
        get_kprocs(6)=-1
        get_kprocs(7)=-1

    end function




!!plane wave functions!!!!!!
    function new_g(k,vec)
        type (pwset)                                                        :: new_g

        real(dp), dimension(3),                         intent(in)          :: k
        integer, dimension(:,:),                        intent(in)          :: vec

        integer                                                             :: ipw

        new_g%k=k

        ALLOCATE( new_g%vec(3,size(vec,2)), &
                  new_g%ind(minval(vec(1,:)):maxval(vec(1,:)), &
                            minval(vec(2,:)):maxval(vec(2,:)), &
                            minval(vec(3,:)):maxval(vec(3,:))) )
        new_g%vec=vec

        new_g%ind=0
        forall (ipw=1:size(vec,2))    new_g%ind(vec(1,ipw),vec(2,ipw),vec(3,ipw))=ipw
    end function

    function copy_g(g)
        type(pwset)                                                          :: copy_g
        type(pwset),intent(in)                                               :: g
       
        copy_g%k=g%k
 
        ALLOCATE( copy_g%vec(3,size(g%vec,2)), &
                  copy_g%ind(lbound(g%ind,1):ubound(g%ind,1), &
                             lbound(g%ind,2):ubound(g%ind,2), &
                             lbound(g%ind,3):ubound(g%ind,3)) )
        copy_g%vec=g%vec
        copy_g%ind=g%ind
 
    end function


    function new_g_ref(k,vec)
        type (pwset),pointer                                                :: new_g_ref

        real(dp), dimension(3),                         intent(in)          :: k
        integer, dimension(:,:),                        intent(in)          :: vec

        integer                                                             :: ipw


        allocate(new_g_ref)
        ALLOCATE( new_g_ref%vec(3,size(vec,2)), &
                  new_g_ref%ind(minval(vec(1,:)):maxval(vec(1,:)), &
                                minval(vec(2,:)):maxval(vec(2,:)), &
                                minval(vec(3,:)):maxval(vec(3,:))) )

        new_g_ref%k=k 
        new_g_ref%vec=vec
        new_g_ref%ind=0
        forall (ipw=1:size(vec,2))    new_g_ref%ind(vec(1,ipw),vec(2,ipw),vec(3,ipw))=ipw
    end function
 

    function copy_g_ref(g)
        type(pwset), pointer                                                 :: copy_g_ref
        type(pwset),intent(in)                                               :: g


        
        allocate(copy_g_ref)
        ALLOCATE( copy_g_ref%vec(3,size(g%vec,2)), &
                  copy_g_ref%ind(lbound(g%ind,1):ubound(g%ind,1), &
                                 lbound(g%ind,2):ubound(g%ind,2), &
                                 lbound(g%ind,3):ubound(g%ind,3)) )

        copy_g_ref%k=g%k
        copy_g_ref%vec=g%vec
        copy_g_ref%ind=g%ind
 
    end function



 

!!type wf functions!!!!!!

    subroutine destroy_wf(wf)
        type(wavefunction)                                                  :: wf

        if (allocated(wf%coeff)) deallocate(wf%coeff)
        if (associated(wf%g))  deallocate(wf%g)

    end subroutine


    function new_wf_vec(coeff, k,gvec)
        type(wavefunction)                                                  :: new_wf_vec

        complex(dpc), dimension(:,:),                   intent(in)          :: coeff
        real(dp), dimension(3),                         intent(in)          :: k
        integer, dimension(:,:),                        intent(in)          :: gvec

        allocate(new_wf_vec%coeff(size(coeff,1),size(coeff,2)))
        new_wf_vec%coeff=coeff

        new_wf_vec%g=>g_(k,gvec)
        
    end function

    function new_wf_ref(coeff, g)
        type(wavefunction)                                                  :: new_wf_ref

        complex(dpc), dimension(:,:),                   intent(in)          :: coeff
        type(pwset), pointer,                           intent(in)          :: g

!        if (size(g%vec,2) /= size(coeff,1)) then
!            write(6,*) "g,coeff mismatch"
!        endif
        
        allocate(new_wf_ref%coeff(size(coeff,1),size(coeff,2)))
        new_wf_ref%coeff=coeff

        new_wf_ref%g=>g

    end function


    function copy_wf_ref(wf)
        type(wavefunction)                                                  :: copy_wf_ref
        type(wavefunction)                                                  :: wf


!        if (size(g%vec,2) /= size(coeff,1)) then
!            write(6,*) "g,coeff mismatch"
!        endif
        
        allocate(copy_wf_ref%coeff(size(wf%coeff,1),size(wf%coeff,2)))
        copy_wf_ref%coeff=wf%coeff

        copy_wf_ref%g=>wf%g

    end function

    function sym_op(wf,ksym,lsym,spinsym)
      USE constants, ONLY : ci, pi
        type(wavefunction)                                                  :: sym_op
        type(wavefunction),                            intent(in)           :: wf
        integer, dimension(3,3),                       intent(in)           :: ksym
        real(dp), dimension(3),                        intent(in)           :: lsym
        complex(dpc), dimension(2,2),                  intent(in), optional :: spinsym

        integer                                                             :: i

        sym_op=reorder_basis(new_wf(wf%coeff, g_(matmul(ksym,wf%g%k),matmul(transpose(ksym),wf%g%vec))),wf%g)
        FORALL (i = 1:size(sym_op%g%vec,2)) &
            sym_op%coeff(i,:) = sym_op%coeff(i,:) * exp(ci*2*pi * dot_product(lsym,sym_op%g%vec(:,i)))

        if (present(spinsym).and.size(sym_op%coeff,2)==2) sym_op%coeff=matmul(sym_op%coeff,transpose(spinsym))
    end function

    function bz_wrap(wf, kshift)
        type(wavefunction)                                                  :: bz_wrap
        type(wavefunction),                             intent(in)          :: wf
        integer, dimension(3),                          intent(in)          :: kshift


        integer, dimension(3,size(wf%g%vec,2))                              :: vec
        integer                                                             :: i

        forall (i=1:size(vec,2)) vec(:,i)=wf%g%vec(:,i)-kshift
        
        bz_wrap=new_wf(wf%coeff,g_(wf%g%k+kshift,vec))
!        bz_wrap=reorder_basis(new_wf(wf%coeff,g_(vec)),wf%g)

    end function

    
  
    pure logical function safe_wv(wf,wv)
        type(wavefunction),                            intent(in)           :: wf
        integer,                     dimension(3),     intent(in)           :: wv

        safe_wv=.false.
        if (all(wv <= ubound(wf%g%ind)) .and. all(wv >= lbound(wf%g%ind)))  safe_wv=(wf%g%ind(wv(1),wv(2),wv(3))>0)  
    end function


    function reorder_basis(wf,g)
        type(wavefunction)                                                  :: reorder_basis

        type(wavefunction),                            intent(in)           :: wf
        type(pwset),                                   intent(in)           :: g


        integer                                                             :: ipw

        allocate(reorder_basis%coeff(size(g%vec,2),size(wf%coeff,2)))

        reorder_basis%coeff=0d0
        FORALL(ipw = 1:size(g%vec,2),safe_wv(wf,g%vec(:,ipw))) &
            reorder_basis%coeff(ipw,:) = wf%coeff(wf%g%ind(g%vec(1,ipw),g%vec(2,ipw),g%vec(3,ipw)),:)

        reorder_basis%g=>g_(g)
        reorder_basis%g%k=wf%g%k
               
    end function


    pure  function inner_prod(wf1,wf2)
        type(wavefunction),                            intent(in)           :: wf1,wf2
        complex(dpc)                                                        :: inner_prod
                
        inner_prod=SUM(CONJG(wf1%coeff)*wf2%coeff)

    end function


    pure function p_overlap(wf1,wf2)
    ! transition dipole matrix elements
        complex(dpc), dimension(3)                                          :: p_overlap
        type(wavefunction),                            intent(in)           :: wf1,wf2
        


        p_overlap=matmul(wf2%g%vec,sum(conjg(wf1%coeff)*wf2%coeff,2))

    end function 
    
!    pure function momentum_overlap(wf1,wf2,kpt)
!    end function
       
    integer function band_block(kpt, ibnd) 
        type(kpt_data),                                 intent(in)          :: kpt
        integer,                                        intent(in)          :: ibnd

        band_block=1
        do while(ibnd > kpt%block(band_block)%band_range(2))
            band_block=band_block+1
        enddo

    end function

    subroutine ft_wf(rgrid,wf,gridsize)
        complex(dpc), dimension(:,:,:,:), allocatable, intent(out)          :: rgrid
        
        type(wavefunction),                             intent(in)          :: wf
        integer, dimension(3),                          intent(in)          :: gridsize


        integer                                                             :: ipw, ispin
        real(dp), dimension(size(wf%coeff,2))                               :: phase

        integer(kind=8)                                                     :: plan

        
        allocate(rgrid(gridsize(1),gridsize(2),gridsize(3),size(wf%coeff,2)))


        rgrid=0d0
        do ispin=1,size(wf%coeff,2)
            do ipw=1,size(wf%g%vec,2)
                rgrid( modulo(wf%g%vec(1,ipw),gridsize(1))+1, &
                       modulo(wf%g%vec(2,ipw),gridsize(2))+1, &
                       modulo(wf%g%vec(3,ipw),gridsize(3))+1, &
                       ispin) = wf%coeff(ipw,ispin)
            enddo
        enddo
       
 
        phase=trwf_phase(wf)
        do ispin=1,size(wf%coeff,2)
            rgrid(:,:,:,ispin)=exp(dcmplx(0d0,-phase(ispin)))*rgrid(:,:,:,ispin)
        enddo
        
        


        
        do ispin=1,size(wf%coeff,2)
            call dfftw_plan_dft_3d(plan,gridsize(1),gridsize(2),gridsize(3),rgrid(:,:,:,ispin),rgrid(:,:,:,ispin),-1,64)
            call dfftw_execute_dft(plan, rgrid(:,:,:,ispin), rgrid(:,:,:,ispin))
            call dfftw_destroy_plan(plan)
        enddo
    
    end subroutine

    function trwf_phase(wf)
        type(wavefunction),                             intent(in)          :: wf

        real(dp), dimension(size(wf%coeff,2))                               :: trwf_phase
        integer                                                             :: ipw, ispin,nlarge

        real, parameter                                                     :: tol=1d-4

        trwf_phase=0d0

        if (any(abs(wf%g%k)>tol)) return 
        
        do ispin=1,size(wf%coeff,2)
            nlarge=count(abs(wf%coeff(:,ispin))>tol)
            do ipw=1,size(wf%coeff,1)
                if (abs(wf%coeff(ipw,ispin))>tol)  trwf_phase=trwf_phase+modulo(aimag(log(wf%coeff(ipw,ispin))),pi)/nlarge
            enddo
        enddo
    end function





    subroutine print_wf(uid, wf,wf2)
        integer,                                        intent(in)          :: uid
        type(wavefunction),                             intent(in)          :: wf
        type(wavefunction),                             intent(in),optional :: wf2
        
        type(wavefunction)                                                  :: wf_reorder
        integer                                                             :: ipw

        if (present(wf2)) then
            wf_reorder=reorder_basis(wf2, wf%g)
            write(uid,*) 'k1= ', wf%g%k, 'k2=', wf2%g%k
            write(uid,*) '  nx   ny   nz     Re(u)       Im(u)        Re(u)       Im(u)'
            IF (size(wf%coeff,2) == 1) &
                WRITE(uid, '(3I5, 4F12.6)') (wf%g%vec(:,ipw), wf%coeff( ipw,:),wf_reorder%coeff(ipw,:), ipw=1,size(wf%g%vec,2))

            if (size(wf%coeff,2)==2)  then
                do ipw=1,size(wf%g%vec,2)
                    WRITE(uid, '(3I5, 6F12.6)') wf%g%vec(:,ipw), wf%coeff(ipw,1), &
                        wf_reorder%coeff(ipw,1), wf%coeff(ipw,1)/wf_reorder%coeff(ipw,1)
                    WRITE(uid, '(A, 6F12.6)') '               ', wf%coeff(ipw,2), &
                        wf_reorder%coeff(ipw,2), wf%coeff(ipw,2)/wf_reorder%coeff(ipw,2)
                enddo
            endif

        else
            write(uid,*) 'k1= ', wf%g%k
            write(uid,*) '  nx   ny   nz     Re(u)       Im(u)'
            if (size(wf%coeff,2)==1)        write(uid,'(3I5, 2F12.6)')(wf%g%vec(:,ipw), wf%coeff( ipw,:),ipw=1,size(wf%g%vec,2))

            if (size(wf%coeff,2)==2)  then
                do ipw=1,size(wf%g%vec,2)
                    write(uid,'(3I5,2F12.6)')wf%g%vec(:,ipw), wf%coeff( ipw,1)
                    write(uid,'(A, 2F12.6)')'               ',wf%coeff( ipw,2)
                enddo
            endif
        endif


    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
