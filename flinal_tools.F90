!!could use blas/lapack of course, but this allows purely functional procedures to be designed, enhancing readability and simplifying debugging.  also,
!!pedagogically edifying for me to actually do

module flinal_tools

  USE constants, ONLY : dp, dpc, eq_tol
implicit none


interface trace
    module procedure ctrace,rtrace
end interface

interface det
    module procedure rdet,cdet
end interface

interface lup
    module procedure rlup,clup
end interface

interface inverse
    module procedure rinverse
end interface

interface lower_inverse
    module procedure rlinverse
end interface

interface cholesky
    module procedure rcholesky
end interface


contains
!!basic procedures!!!!!!!!!!!!!!!!!!!!!!
    elemental complex(dpc) function safe_clog(x)
        complex(dpc),                                           intent(in)          :: x
    
        if (dble(x)>eq_tol) then
            if (abs(dimag(x)/dble(x))<10.0) then
            
                safe_clog=log(x)
            endif
        else
            safe_clog=0d0
        endif
    end function

    pure function tranjg(m)
        complex(dpc), dimension(:,:),                               intent(in)      :: m
        complex(dpc), dimension(size(m,2),size(m,1))                                :: tranjg
        

        tranjg=conjg(transpose(m))
    end function

    function mat_log(matrix)
        complex(dpc), dimension(:,:),                               intent(in)      :: matrix

        complex(dpc), dimension(size(matrix,1),size(matrix,2))                      :: mat_log


        complex(dpc), dimension(:,:),                   allocatable                 :: reig_vec,leig_vec,d
        integer                                                                     :: i,ndim

        ndim=size(matrix,1)

        allocate(reig_vec(ndim,ndim),leig_vec(ndim,ndim),d(ndim,ndim))
    

        call diagonalize(matrix,d,reig_vec,leig_vec)

        mat_log=matmul(matmul(reig_vec,safe_clog(d)),tranjg(leig_vec))
    end function

    subroutine diagonalize(matrix,d,runitary,lunitary)
        complex(dpc), dimension(:,:),                               intent(in)      :: matrix
        complex(dpc), dimension(:,:),                               intent(out)     :: d, runitary, lunitary

        complex(dpc), dimension(:,:),                   allocatable                 :: reig_vec,leig_vec, temp_diag
        complex(dpc), dimension(:),                     allocatable                 :: eig_val

        integer                                                                     :: i,ndim

        ndim=size(d,1)

        allocate(reig_vec(ndim,ndim),leig_vec(ndim,ndim),eig_val(ndim),temp_diag(ndim,ndim))

        call eigen(matrix,eig_val,reig_vec,leig_vec)
        
        temp_diag=MATMUL(tranjg(leig_vec),reig_vec)

        d=0d0
        forall (i=1:ndim)
            d(i,i)=eig_val(i)
            runitary(:,i)=reig_vec(:,i)/sqrt(temp_diag(i,i))
            lunitary(:,i)=leig_vec(:,i)/conjg(sqrt(temp_diag(i,i)))
        end forall
        

    end subroutine


    complex(dpc) function zdet(m)
        complex(dpc), DIMENSION(:,:),               intent(in)  :: m

        complex(dpc), dimension(size(m,1),size(m,2))     :: m_temp
        INTEGER                                          :: i,ndim

        INTEGER                                          :: lwork,info
        INTEGER, ALLOCATABLE, DIMENSION(:)               :: ipiv
        complex(dpc), ALLOCATABLE, DIMENSION(:)          :: work

        IF (size(m,1).ne.size(m,2)) then
            zdet=0d0
            RETURN
        ELSE
            ndim=size(m,1)
            lwork=4*ndim

            ALLOCATE(work(lwork),ipiv(ndim))
            m_temp=m
            call zgetrf(ndim,ndim,m_temp,ndim,ipiv,info)

            zdet=1d0
            do i=1,ndim
                zdet=zdet*m_temp(i,i)
            enddo

            DEALLOCATE(work,ipiv)
        ENDIF
            
    end function

!    FUNCTION inverse(m)
!        complex(dpc), DIMENSION(:,:)                     :: m
!        complex(dpc), DIMENSION(size(m,1),size(m,1))     :: inverse
!        INTEGER                                          :: ndim
!
!        INTEGER                                          :: lwork,info
!        INTEGER, ALLOCATABLE, DIMENSION(:)               :: ipiv
!        complex(dpc), ALLOCATABLE, DIMENSION(:)          :: work
!
!
!
!        IF (size(m,1).ne.size(m,2)) then
!            inverse=0d0
!            RETURN
!        ELSE
!            ndim=size(m,1)
!            lwork=4*ndim
!
!            ALLOCATE(work(lwork),ipiv(ndim))
!            inverse=m
!            call zgetrf(ndim,ndim,inverse,ndim,ipiv,info)
!            call zgetri(ndim,inverse,ndim,ipiv,work,lwork, info)

!            DEALLOCATE(work,ipiv)
!        ENDIF
!    END FUNCTION inverse

        
    subroutine eigen(m,eigval,reigvec,leigvec)
        complex(dpc), dimension(:,:),                           intent(in)    :: m
        complex(dpc), dimension(:,:),                           intent(out)   :: reigvec
        complex(dpc), dimension(:),                             intent(out)   :: eigval
        complex(dpc), dimension(:,:),               optional,   intent(out)   :: leigvec
      
        !!!Indices
        integer                                                               :: i,j


        !!!ZGEEV
        integer                                                               :: ipiv(size(eigval)),info,lwork
        complex(dpc), dimension(size(eigval,1))                               :: w
        complex(dpc), dimension(size(reigvec,1),size(reigvec,2))              :: m_temp,vl,vr
        complex(dpc), dimension(6*size(eigval))                               :: work
        real(dp),     dimension(2*size(eigval))                               :: rwork
        
  
        m_temp=m
        !!!Solve eigenvalue equation                              
        call zgeev('V','V',size(eigval),m_temp,size(eigval),w,vl,size(eigval),vr,size(eigval),work,size(work),rwork,info)

        reigvec=vr                                                            
        eigval=w
        if (present(leigvec)) leigvec=vl                                      

    end subroutine

    pure real(dp) function rtrace(matrix)
        real(dp), dimension(:,:),                               intent(in)      :: matrix

        integer                                                                 :: i

        rtrace=0d0
        do i=1,size(matrix,1)
            rtrace=rtrace+matrix(i,i)
        enddo

    end function

    pure complex(dpc) function ctrace(matrix)
        complex(dpc), dimension(:,:),                           intent(in)      :: matrix

        integer                                                                 :: i

        ctrace=0d0
        do i=1,size(matrix,1)
            ctrace=ctrace+matrix(i,i)
        enddo

    end function




    !!cholesky decomposes a real (posdef) matrix.  returns the lower triangular
    pure function rcholesky(m)
        real(dp),dimension(:,:),                                intent(in)      :: m
        
        real(dp),dimension(size(m,1),size(m,2))                                 :: rcholesky

        integer                                                                 :: i,j

        rcholesky=0d0

        !!first column
        rcholesky(1,:)=m(1,:)/sqrt(m(1,1))
        
        !!2-n column
        do i=2,size(m,1)
            do j=i,size(m,2)
                rcholesky(i,j)=m(i,j)-dot_product(rcholesky(1:i-1,i),rcholesky(1:i-1,j))
            enddo
            rcholesky(i,:)=rcholesky(i,:)/sqrt(rcholesky(i,i))
        enddo        

        rcholesky=transpose(rcholesky)
        
    end function


    !!does an inverse based on lup decomposition
    pure function rinverse(m)
        real(dp), dimension(:,:),                                   intent(in)      :: m
    
        real(dp), dimension(size(m,1),size(m,2))                                    :: rinverse
    
    
    
        real(dp), dimension(size(m,1),3*size(m,2))                                  :: lup
        integer                                                                     :: n
    
    
        n=size(m,1)
    
        lup=rlup(m)
    
        !!m^-1=P^T,U^-1,L^-1
        rinverse=matmul(transpose(lup(:,2*n+1:3*n)),matmul(transpose(rlinverse(transpose(lup(:,n+1:2*n)))),rlinverse(lup(:,:n))))
    
    
    end function
    
    !! inverts a lower triangular matrix via back substitution
    pure function rlinverse(l)
        real(dp), dimension(:,:),                                   intent(in)      :: l
    
        real(dp), dimension(size(l,1),size(l,2))                                    :: rlinverse
    
    
        integer                                                                     :: i,j,n
    
        n=size(l,1)
        rlinverse=0.0
    
        do i=n,1,-1
            rlinverse(i,i)=1/l(i,i)
            forall (j=i+1:n) rlinverse(j,i)=-dot_product(rlinverse(j,i+1:j),l(i+1:j,i))/l(i,i)
        enddo
    
    end function
    
    
    pure function iidentity(n)
        integer,                                                intent(in)          :: n
    
        integer, dimension(n,n)                                                     :: iidentity
    
    
        integer                                                                     :: i
    
        iidentity=0
        forall (i=1:n) iidentity(i,i)=1
    end function
    
    pure function ridentity(n)
        integer,                                                intent(in)          :: n
    
        real(dp), dimension(n,n)                                                    :: ridentity
    
    
        integer                                                                     :: i
    
        ridentity=0d0
        forall (i=1:n) ridentity(i,i)=1d0
    end function
    
    pure FUNCTION cidentity(n)
        integer,                                                intent(in)          :: n
    
        complex(dpc), dimension(n,n)                                                :: cidentity
    
    
        integer                                                                     :: i
    
        cidentity=(0d0,0d0)
        forall (i=1:n) cidentity(i,i)=(1d0,0d0)
    end function
 
    !!lup decomposition
    pure function clup(a)
        complex(dpc), dimension(:,:),                               intent(in)      :: a
    
        complex(dpc), dimension(size(a,1),3*size(a,2))                              :: clup
    
    
        complex(dpc), dimension(size(a,1))                                          :: swap
        integer                                                                     :: i,j,n
    
        n=size(a,1)
    
        clup(:,n+1:2*n)=a
        clup(:,:n)=ridentity(n)
        clup(:,2*n+1:)=ridentity(n)
    
        do i=1, n
            if (abs(clup(i,i+n)) <= eq_tol) then
                j=1
                do while(abs(clup(i,i+n+j)) <= eq_tol)
                    j=j+1
                enddo
                swap=clup(:,i+n)
                clup(:,i+n)=clup(:,i+n+j)
                clup(:,i+n+j)=swap
    
                swap=clup(:,i+2*n)
                clup(:,i+2*n)=clup(:,i+2*n+j)
                clup(:,i+2*n+j)=swap
            endif
    
            clup(i+1:,i)=-clup(i+1:,i+n)/clup(i,i+n)
            clup(i:,i+n:2*n)=matmul(clup(i:,i:n),clup(i:,i+n:2*n))
        enddo
    
        clup(:,:n)=-clup(:,:n)+2*ridentity(n)
    end function
    

    pure function rlup(a)
        real(dp), dimension(:,:),                                   intent(in)      :: a
    
        real(dp), dimension(size(a,1),3*size(a,2))                                  :: rlup
    
    
        real(dp), dimension(size(a,1))                                              :: swap
        integer                                                                     :: i,j,n
    
        n=size(a,1)
    
        rlup(:,n+1:2*n)=a
        rlup(:,:n)=ridentity(n)
        rlup(:,2*n+1:)=ridentity(n)
    
        do i=1, n
            if (abs(rlup(i,i+n)) <= eq_tol) then
                j=1
                do while(abs(rlup(i,i+n+j)) <= eq_tol)
                    j=j+1
                enddo
                swap=rlup(:,i+n)
                rlup(:,i+n)=rlup(:,i+n+j)
                rlup(:,i+n+j)=swap
    
                swap=rlup(:,i+2*n)
                rlup(:,i+2*n)=rlup(:,i+2*n+j)
                rlup(:,i+2*n+j)=swap
            endif
    
            rlup(i+1:,i)=-rlup(i+1:,i+n)/rlup(i,i+n)
            rlup(i:,i+n:2*n)=matmul(rlup(i:,i:n),rlup(i:,i+n:2*n))
        enddo
    
        rlup(:,:n)=-rlup(:,:n)+2*ridentity(n)
    end function
    

    pure complex(dpc) function cdet(m)
        complex(dpc),dimension(:,:),                            intent(in)      :: m

        complex(dpc), dimension(size(m,1),3*size(m,2))                          :: l
        complex(dpc), dimension(size(m,1))                                      :: temp
        integer                                                                 :: i,j,n

        n=size(m,1)
        
        l=clup(m)
  
        cdet=1d0
        do i=1,n
            j=2*n+i
            do while(l(i,j)/=1.0)
                j=j+1
            enddo
            if (i/=j) then
                temp=l(:,i+2*n)
                l(:,i+2*n)=l(:,j)
                l(:,j)=temp
                cdet=-cdet
            endif
        enddo
        
            
        do i=1,n
            cdet=cdet*l(i,n+i)
        enddo


    end function            

    !!Get's sign wrong!! ok, for now, since only use needs just abs val, but must be fixed
    pure real(dp) function rdet(m)
        real(dp),dimension(:,:),                                intent(in)      :: m

        real(dp), dimension(size(m,1),3*size(m,2))                              :: l
        real(dp), dimension(size(m,1))                                          :: temp
        integer                                                                 :: i,j,n

        n=size(m,1)
        
        l=rlup(m)
  
        rdet=1d0
        do i=1,n
            j=2*n+i
            do while(l(i,j)/=1.0)
                j=j+1
            enddo
            if (i/=j) then
                temp=l(:,i+2*n)
                l(:,i+2*n)=l(:,j)
                l(:,j)=temp
                rdet=-rdet
            endif
        enddo
        
            
        do i=1,n
            rdet=rdet*l(i,n+i)
        enddo


    end function            

!    real(dp) function det(m)
!        real(dp),dimension(:,:),                                intent(in)      :: m
!
!        real(dp), dimension(size(m,1),size(m,2))                                :: l
!        integer                                                                 :: i
!        
!        l=cholesky(m)
!  
!        det=1d0
!        do i=1,size(m,1)
!            det=det*l(i,i)**2
!        enddo
!            
!    end function

!    function rinverse(m)
!        real(dp),dimension(:,:),                                intent(in)      :: m
!        
!        real(dp),dimension(size(m,1),size(m,2))                                 :: L,rinverse
!
!        
!        integer                                                                 :: i,j

!        L=cholesky(m)
!        write(6,*)L

!        rinverse=0d0
!        do i=1,size(m,2)
!            rinverse(i,i)=1/L(i,i)
!            do j=i+1,size(m,1)
!                rinverse(j,i)=-dot_product(L(j,1:j-1),rinverse(1:j-1,i))/L(j,j)
!            enddo
!        enddo

!        rinverse=matmul(transpose(rinverse),rinverse)


!    end function 


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
