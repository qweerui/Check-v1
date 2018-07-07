module Mathmatic
!!===============================================================================84
!!
!!  Simplified edition of MKL/Lapack lib
!!      
!!  future work :
!!      1.squareMarix().
!!      2.Unified _r4 and _r8 verson.
!!
!!      contains:
!!              sub:
!!
!!                  SVD_Lapack_r4(m,n,AA,UU,SigMa,VVT)
!!                  SVD_Lapack_r8(m,n,AA,UU,SigMa,VVT)
!!                  SVD_BasedOnEOF_r8(m,n,AA,UU,SigMa,VVT)
!!
!!                                      **** SVD probrems:
!!                                      UU /= VVT  UUi not always eq VVTi 
!!                                      should be used with extreme caution.
!!                                      **** above problem solved 
!!                                      AA is not Symmetric. AA = AAT
!!                  !-----------------------------------------------------------!
!!                  !SVD:                                                       !
!!                  !     A   =  U  * SigMa *  VT                               !
!!                  !    mxn    mxm    mxn    nxn                               !
!!                  !                                                           !
!!                  !SVD_Lapack:                                                !
!!                  !     AA  = UU  * SigMa *  VVT                              !
!!                  !    mxn    mxn    nxn    nxn                               !
!!                  !                                                           !
!!                  !SVD_BasedOnEOF:                                            !
!!                  !     BB  = AAT * AA                                        !
!!                  !     BB  = VV * SigMaE * VVT   (EOF)                       !
!!                  !                                                           !
!!                  !     AA  = AA * I = AA * VV*VVT                            !
!!                  !         = {AA*VV*sqrt(SigMaE)-1} * sqrt(SigMaE) * VVT     !
!!                  !    mxn             mxn               nxn          nxn     !
!!                  !     AA  =          UU            *  SigMa       * VVT     !
!!                  !                                                           !
!!                  !-----------------------------------------------------------!
!!
!!                  EOF_Lapack_r4(n,AA,UU,SigMa)
!!                  EOF_Lapack_r8(n,AA,UU,SigMa)
!!
!!                                           ONLY:
!!                                      **** Symmetric Eigenvalue Problems !!!
!!                                      **** Checked!
!!
!!                  PRINT_MATRIX_r4(name,m,A)
!!                  PRINT_MATRIX_r8(name,m,A)
!!
!!                                      **** Checked!
!!
!!                  QRsquare_r4(m,AA,QQ,RR)
!!                  QRsquare_r8(m,AA,QQ,RR)
!!              
!!                                      **** Checked!
!!
!!              function:
!!                  inv_r8(A) result(Ainv)
!!
!!                              **** can be precise even n = 2000 or greater
!!
!!                  inv_r4s(A) result(Ainv)
!!                  inv_r4(A) result(Ainv)
!!
!!                              **** while using matmul when n > 200, 
!!                              **** truncation error can not be ignored.
!!
!!                  matmul(A,A)
!!                      c(m,l) = matmul(a(m,n),b(n,l))
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2016.05.25               Hanrui           original code
!!
!!        2016.06.21               Hanrui           QR: QRsquare_r8
!!                                                      QRsquare_r4
!!
!!        2016.06.23               Hanrui           EOF: Symmetric only. checked
!!                                  see: /home/ruihan/example/MKL_test/001/QR.f90
!!
!!        2016.06.23               Hanrui           SVD: checked
!!                                  see: /home/ruihan/example/MKL_test/001/QR.f90
!!
!!        2016.06.26               Hanrui           SVD: SVD_BasedOnEOF_r8
!!                                                  use only m >> n
!!
!!===============================================================================84

    implicit none
    !--------------------------------------------------------------------------
    !   inverse matrix
    PRIVATE     ::  inv_r4s
    PRIVATE     ::  inv_r4
    !
    PUBLIC      ::  inv_r8

    !   SVD decomposition
    PRIVATE     ::  SVD_Lapack_r4
    !
    PUBLIC      ::  SVD_Lapack_r8
    PUBLIC      ::  SVD_BasedOnEOF_r8

    !   EOF decomposition
    PRIVATE     ::  EOF_Lapack_r4
    !
    PUBLIC      ::  EOF_Lapack_r8

    !   QR decomposition
    PRIVATE     ::  QRsquare_r4
    !
    PUBLIC      ::  QRsquare_r8
    
    !   print matrix
    PRIVATE     ::  PRINT_MATRIX_r4
    !
    PUBLIC      ::  PRINT_MATRIX_r8

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function inv_r4s(A) result(Ainv)
    !!=======================================================================76
    !! Returns the inverse of a matrix calculated by finding the LU
    !! decomposition.  Depends on LAPACK.
    !!
    !!  with greater error to n>200
    !!  comes from truncation error while using matmul
    !!=======================================================================76
        
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real, dimension(:,:), intent(in)       :: A
            real, dimension(size(A,1),size(A,2))   :: Ainv

            real(8), dimension(size(A,1),size(A,2))   :: A_r8
            real(8), dimension(size(A,1),size(A,2))   :: Ainv_r8

            !----------------------------------------------------------
            !   parameter
            real(8), dimension(size(A_r8,1)) :: work  ! work array for LAPACK
            integer, dimension(size(A_r8,1)) :: ipiv   ! pivot indices
            integer :: n, info

            !----------------------------------------------------------

            ! External procedures defined in LAPACK
            external DGETRF
            external DGETRI
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------

            A_r8 = real( A  , kind = 8 )

            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv_r8 = A_r8
            n = size(A_r8,1)

            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call DGETRF(n, n, Ainv_r8, n, ipiv, info)

            if (info /= 0) then
               stop 'Matrix is numerically singular!'
            if(info.lt.0)then
                print '(" lu decomposition:  illegal value ")'
                stop
            endif
            if(info.gt.0)then
                write(*,"('lu decomposition: u(',i4,',',i4,') = 0 ')")info,info
            endif
            end if

            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call DGETRI(n, Ainv_r8, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed!'
            end if

            Ainv_r8 = matmul(A_r8,Ainv_r8)
            Ainv = real( Ainv_r8  , kind = 4 )
            !----------------------------------------------------------

        !------------------------------------------------------------------
            
    !!=======================================================================76
    end function inv_r4s

    function inv_r4(A) result(Ainv)
    !!=======================================================================76
    !! Returns the inverse of a matrix calculated by finding the LU
    !! decomposition.  Depends on LAPACK.
    !!  
    !!  with greater error to n>200
    !!  comes from truncation error while using matmul
    !!=======================================================================76
        
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real, dimension(:,:), intent(in)       :: A
            real, dimension(size(A,1),size(A,2))   :: Ainv

            !----------------------------------------------------------
            !   parameter
            real, dimension(size(A,1))    :: work  ! work array for LAPACK
            integer, dimension(size(A,1)) :: ipiv   ! pivot indices
            integer :: n, info
            !----------------------------------------------------------

            ! External procedures defined in LAPACK
            external SGETRF
            external SGETRI
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------

            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv = A
            n = size(A,1)

            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call SGETRF(n, n, Ainv, n, ipiv, info)

            if (info /= 0) then
               stop 'Matrix is numerically singular!'
            if(info.lt.0)then
                print '(" lu decomposition:  illegal value ")'
                stop
            endif
            if(info.gt.0)then
                write(*,"('lu decomposition: u(',i4,',',i4,') = 0 ')")info,info
            endif
            end if

            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call SGETRI(n, Ainv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed!'
            end if
            !----------------------------------------------------------

        !------------------------------------------------------------------
            
    !!=======================================================================76
    end function inv_r4

    function inv_r8(A) result(Ainv)
    !!=======================================================================76
    !! Returns the inverse of a matrix calculated by finding the LU
    !! decomposition.  Depends on LAPACK.
    !!=======================================================================76
        
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real(8), dimension(:,:), intent(in)       :: A
            real(8), dimension(size(A,1),size(A,2))   :: Ainv

            !----------------------------------------------------------
            !   parameter
            real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
            integer, dimension(size(A,1)) :: ipiv   ! pivot indices
            integer :: n, info
            !----------------------------------------------------------

            ! External procedures defined in LAPACK
            external DGETRF
            external DGETRI
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------

            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv = A
            n = size(A,1)

            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call DGETRF(n, n, Ainv, n, ipiv, info)

            if (info /= 0) then
               stop 'Matrix is numerically singular!'
            if(info.lt.0)then
                print '(" lu decomposition:  illegal value ")'
                stop
            endif
            if(info.gt.0)then
                write(*,"('lu decomposition: u(',i4,',',i4,') = 0 ')")info,info
            endif
            end if

            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call DGETRI(n, Ainv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed!'
            end if
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end function inv_r8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine SVD_Lapack_r4(m,n,AA,UU,SigMa,VVT)
    !!=======================================================================76
    !!give a single value decomposition of matrix A:
    !!  A   =  U  * SigMa *  VT 
    !! mxn    mxm    mxn    nxn
    !!
    !! outputs:
    !!  AA  = UU  * SigMa *  VVT
    !! mxn    mxn    nxn    nxn
    !!
    !!use MKL Lapack lib
    !! contains: 
    !!      dgbbrd    P704  MKL man
    !!      dbdsqr    P722  MKL man
    !!      dgesvd   P1093  MKL man
    !!
    !!  tmp = matmul(UU,SigMa)
    !!  tmp = matmul(tmp,VVT) - Abak
    !!  call PRINT_MATRIX_r8("UU*sigma*VVT - A",m,tmp)  tmp = 0.0~0.0
    !!
    !!  UU(:,i) .O. UU(:,i)
    !!  Left singular vectors (stored columnwise)
    !!  VVT(i,:) .O. VVT(i,:)
    !!  Right singular vectors (stored rowwise)     Checked !!! 
    !!  see: /home/ruihan/example/MKL_test/001/QR.f90
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  m,  n  
            real,intent(in)        ::  AA(m,n)
            real,intent(inout)     ::  UU(m,n)
            real,intent(inout)     ::  SigMa(n,n)
            real,intent(inout)     ::  VVT(n,n)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  ldu
            integer                ::  ldvt
            integer                ::  lwmax
            parameter                  (lwmax = 100000)

            !   local scalars
            integer                ::  info,   lwork

            !   local arrays
            real,allocatable       ::  A(:,:)
            real,allocatable       ::  U(:,:)
            real,allocatable       ::  Vt(:,:)
            real                   ::  S(n)
            real                   ::  work(lwmax)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            lda  = m
            ldu  = m
            ldvt = n

            allocate(A(lda,n))
            allocate(U(ldu,m))
            allocate(Vt(ldvt,m))

            A = AA

            !   Executable Statements
            write(*,*)'dgesvd :the arithmetic operation of SVD'

            !   Query the optimal workspace.
            lwork = -1
            call sgesvd('A', 'A', m, n, A, lda, S, U, ldu, Vt, ldvt, work, lwork, info)
            lwork = min( lwmax, int( work(1) ) )

            !   Compute SVD.
            call sgesvd('A', 'A', m, n, A, lda, S, U, ldu, Vt, ldvt, work, lwork, info)

            !   Check for convergence.
            if( info.gt.0 ) then
            write(*,*) "info: ",info
            write(*,*) "work(1): ",work(1)
            write(*,*) 'the algorithm computing SVD failed to converge.'
            stop
            end if  

            !   assign SigMa
            SigMa = 0.0
            do i=1,n
            SigMa(i,i) = S(i)
            enddo

            UU(:,:)  = U(:,1:n)
            VVT = Vt

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine SVD_Lapack_r4

    subroutine SVD_Lapack_r8(m,n,AA,UU,SigMa,VVT)
    !!=======================================================================76
    !!give a single value decomposition of matrix A:
    !!  A   =  U  * SigMa *  VT 
    !! mxn    mxm    mxn    nxn
    !!
    !! outputs:
    !!  AA  = UU  * SigMa *  VVT
    !! mxn    mxn    nxn    nxn
    !!
    !!use MKL Lapack lib
    !! contains: 
    !!      dgbbrd    P704  MKL man
    !!      dbdsqr    P722  MKL man
    !!      dgesvd   P1093  MKL man
    !!
    !!  tmp = matmul(UU,SigMa)
    !!  tmp = matmul(tmp,VVT) - Abak
    !!  call PRINT_MATRIX_r8("UU*sigma*VVT - A",m,tmp)  tmp = 0.0~0.0
    !!
    !!  UU(:,i) .O. UU(:,i)
    !!  Left singular vectors (stored columnwise)
    !!  VVT(i,:) .O. VVT(i,:)
    !!  Right singular vectors (stored rowwise)     Checked !!! 
    !!  see: /home/ruihan/example/MKL_test/001/QR.f90
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  m,  n  
            real(8),intent(in)     ::  AA(m,n)
            real(8),intent(inout)  ::  UU(m,n)
            real(8),intent(inout)  ::  SigMa(n,n)
            real(8),intent(inout)  ::  VVT(n,n)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  ldu
            integer                ::  ldvt
            integer                ::  lwmax
            parameter                  (lwmax = 100000)

            !   local scalars
            integer                ::  info,   lwork

            !   local arrays
            real(8),allocatable    ::  A(:,:)
            real(8),allocatable    ::  U(:,:)
            real(8),allocatable    ::  Vt(:,:)
            real(8)                ::  S(n)
            real(8)                ::  work(lwmax)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            lda  = m
            ldu  = m
            ldvt = n

            allocate(A(lda,n))
            allocate(U(ldu,m))
            allocate(Vt(ldvt,m))

            A = AA

            !   Executable Statements
            write(*,*)'dgesvd :the arithmetic operation of SVD'

            !   Query the optimal workspace.
            lwork = -1
            call dgesvd('A', 'A', m, n, A, lda, S, U, ldu, Vt, ldvt, work, lwork, info)
            lwork = min( lwmax, int( work(1) ) )

            !   Compute SVD.
            call dgesvd('A', 'A', m, n, A, lda, S, U, ldu, Vt, ldvt, work, lwork, info)

            !   Check for convergence.
            if( info.gt.0 ) then
            write(*,*) "info: ",info
            write(*,*) "work(1): ",work(1)
            write(*,*) 'the algorithm computing SVD failed to converge.'
            stop
            end if  

            !   assign SigMa
            SigMa = 0.0
            do i=1,n
            SigMa(i,i) = S(i)
            enddo

            UU(:,:)  = U(:,1:n)
            VVT = Vt

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine SVD_Lapack_r8

    subroutine SVD_BasedOnEOF_r8(m,n,AA,UU,SigMa,VVT)
    !!=======================================================================76
    !!  AAT * AA = BB
    !!  BB  = VV * SigMaE * VVT   (EOF)
    !!
    !!  AA  = AA * I = AA * VV*VVT
    !!      = {AA*VV*sqrt(SigMaE)-1} * sqrt(SigMaE) * VVT
    !!  mxn           mxn                 nxn         nxn
    !!  AA  =          UU            *  SigMa       * VVT
    !!
    !!  check: UU.O.UU ?
    !!      checked UU(:,i).O.UU(:,i) . 
    !!              UU(:,i) standarded orthogonal vector
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)      ::  m,  n  
            real(8),intent(in)      ::  AA(m,n)
            real(8),intent(inout)   ::  UU(m,n)
            real(8),intent(inout)   ::  SigMa(n,n)
            real(8),intent(inout)   ::  VVT(n,n)

            !----------------------------------------------------------
            !   local arrays
            real(8)                 ::  BB(n,n)
            real(8)                 ::  VV(n,n)
            real(8)                 ::  SigMaE(n,n)
            real(8)                 ::  SigMaEinv(n,n)
            real(8)                 ::  II(n,n)
            real(8)                 ::  IIinv(n,n)
            !   parameter
            integer                 ::  i,  j,  k
            !----------------------------------------------------------
        
        !------------------------------------------------------------------

            !----------------------------------------------------------
            BB = matmul(transpose(AA),AA)

            !!! call EOF_Lapack_r8(n,BB,VV,SigMaE) 
            !!! VVT = transpose(VV)
            !!! same results but ranking from  min to max of eigenvalue.
            call SVD_Lapack_r8(n,n,BB,VV,SigMaE,VVT)

            do i=1,n
            SigMaE(i,i) = sqrt(SigMaE(i,i))
            enddo
            SigMaEinv = inv_r8(SigMaE)
            UU = matmul(AA,VV)
            UU = matmul(UU,SigMaEinv)

            SigMa = SigMaE

            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine SVD_BasedOnEOF_r8


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine EOF_Lapack_r4(n,AA,UU,SigMa)
    !!=======================================================================76
    !!caculate Empirical Orthogonal Function of matrix A:
    !!  A   =  U  * SigMa *  UT
    !! nxn    nxn    nxn    nxn
    !!use MKL Lapack lib
    !! contains: 
    !!      ssyev       P971  MKL man
    !! ONLY:
    !!      Symmetric Eigenvalue Problems   AA = AAT
    !!
    !!  AA = matmul(UU,SigMa)
    !!  AA = matmul(AA,transpose(UU))
    !!  AA = AA - Abak
    !!  call PRINT_MATRIX_r8("U*sigma*UT - A",m,AA)  results: 0.0~0.0
    !!
    !!  UU(:,i) .O. UU(:,i) Checked !!! 
    !!  Eigenvectors (stored columnwise)
    !!  see: /home/ruihan/example/MKL_test/001/QR.f90
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  n  
            real,intent(in)        ::  AA(n,n)
            real,intent(inout)     ::  UU(n,n)
            real,intent(inout)     ::  SigMa(n,n)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  lwmax
            parameter                  (lwmax = 100000)

            !   local scalars
            integer                ::  info,   lwork
            real                   ::  error

            !   local arrays
            real,allocatable       ::  A(:,:)
            real                   ::  W(n)
            real                   ::  work(lwmax)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            lda  = n
            allocate(A(lda,n))

            !check
            A = AA - transpose(AA)
            error = sum(A)
            if(error.ne.0.0) then
                stop "input matrix AA of EOF_Lapack_r4 is not Symmetric."
            endif

            A = AA

            !   Executable Statements
            write(*,*)'DSYEV :the arithmetic operation of EOF'

            !   Query the optimal workspace.
            lwork = -1
            call ssyev( 'vectors', 'upper', n, A, lda, W, work, lwork, info )
            lwork = min( lwmax, int( work(1) ) )

            !   Compute EOF.
            call ssyev( 'vectors', 'upper', n, A, lda, W, work, lwork, info )

            !   Check for convergence.
            if( info.gt.0 ) then
            write(*,*) "info: ",info
            write(*,*) "work(1): ",work(1)
            write(*,*) 'the algorithm computing EOF failed to converge.'
            stop
            end if

            !   assign SigMa
            SigMa = 0.0
            do i=1,n
            SigMa(i,i) = W(i)
            enddo
            
            UU  = A
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine EOF_Lapack_r4

    subroutine EOF_Lapack_r8(n,AA,UU,SigMa)
    !!=======================================================================76
    !!caculate Empirical Orthogonal Function of matrix A:
    !!  A   =  U  * SigMa *  UT
    !! nxn    nxn    nxn    nxn
    !!use MKL Lapack lib
    !! contains: 
    !!      dsyev       P971  MKL man
    !! ONLY:
    !!      Symmetric Eigenvalue Problems    AA = AAT
    !!
    !!  AA = matmul(UU,SigMa)
    !!  AA = matmul(AA,transpose(UU))
    !!  AA = AA - Abak
    !!  call PRINT_MATRIX_r8("U*sigma*UT - A",m,AA)  results: 0.0~0.0
    !!
    !!  UU(:,i) .O. UU(:,i) Checked !!! 
    !!  Eigenvectors (stored columnwise)
    !!  see: /home/ruihan/example/MKL_test/001/QR.f90
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  n  
            real(8),intent(in)     ::  AA(n,n)
            real(8),intent(inout)  ::  UU(n,n)
            real(8),intent(inout)  ::  SigMa(n,n)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  lwmax
            parameter                  (lwmax = 100000)

            !   local scalars
            integer                ::  info,   lwork
            real(8)                ::  error

            !   local arrays
            real(8),allocatable    ::  A(:,:)
            real(8)                ::  W(n)
            real(8)                ::  work(lwmax)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            lda  = n
            allocate(A(lda,n))

            !check
            A = AA - transpose(AA)
            error = sum(A)
            if(error.ne.0.0) then
                stop "input matrix AA of EOF_Lapack_r8 is not Symmetric."
            endif

            A = AA

            !   Executable Statements
            write(*,*)'DSYEV :the arithmetic operation of EOF'

            !   Query the optimal workspace.
            lwork = -1
            call dsyev( 'vectors', 'upper', n, A, lda, W, work, lwork, info )
            lwork = min( lwmax, int( work(1) ) )

            !   Compute EOF.
            call dsyev( 'vectors', 'upper', n, A, lda, W, work, lwork, info )

            !   Check for convergence.
            if( info.gt.0 ) then
            write(*,*) "info: ",info
            write(*,*) "work(1): ",work(1)
            write(*,*) 'the algorithm computing EOF failed to converge.'
            stop
            end if

            !   assign SigMa
            SigMa = 0.0
            do i=1,n
            SigMa(i,i) = W(i)
            enddo
            
            UU = A
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine EOF_Lapack_r8


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine QRsquare_r4(m,AA,QQ,RR)
    !!=======================================================================76
    !! give a QR factorization of square matrix A(m,m):
    !! powered by LAPACK.
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  m 
            real,intent(in)        ::  AA(m,m)
            real,intent(inout)     ::  QQ(m,m)
            real,intent(inout)     ::  RR(m,m)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  lwork

            !   local scalars
            integer                ::  ifail
            integer                ::  info
            real                   ::  err
            real                   ::  eps
            parameter                  (eps = 10E-6)
            
            !   local arrays
            real                   ::  A(m,m)
            real                   ::  tau(m)
            real,allocatable       ::  work(:)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            external sgeqrf, sorgqr

            lda = m
            lwork = 64*m
            allocate(work(lwork))

            A = AA
            !   Compute the QR factorization of A
            call sgeqrf(m,m,A,lda,tau,work,lwork,info)
            !   Form the leading m columns of Q explicitly
            call sorgqr(m,m,m,A,lda,tau,work,lwork,info)

            QQ = A
            RR = matmul(transpose(QQ),AA)
            A  = matmul(QQ,RR) - AA
            err= sum(A)/(m*m*1.0)
            if(err.gt.eps) then
                write(*,*) "Accuracy does not meet the requirements"
                write(*,*) "error:",err,"contrast with eps:",eps
                stop "check point QRsquare_r4 001"
            endif
            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine QRsquare_r4

    subroutine QRsquare_r8(m,AA,QQ,RR)
    !!=======================================================================76
    !! give a QR factorization of square matrix A(m,m):
    !! powered by LAPACK.
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)     ::  m 
            real(8),intent(in)     ::  AA(m,m)
            real(8),intent(inout)  ::  QQ(m,m)
            real(8),intent(inout)  ::  RR(m,m)

            !----------------------------------------------------------
            !   parameter
            integer                ::  i,  j,  k
            integer                ::  lda
            integer                ::  lwork

            !   local scalars
            integer                ::  ifail
            integer                ::  info
            real(8)                ::  err
            real(8)                ::  eps
            parameter                  (eps = 10E-12)
            
            !   local arrays
            real(8)                ::  A(m,m)
            real(8)                ::  tau(m)
            real(8),allocatable    ::  work(:)
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            external dgeqrf, dorgqr

            lda = m
            lwork = 64*m
            allocate(work(lwork))

            A = AA
            !   Compute the QR factorization of A
            call dgeqrf(m,m,A,lda,tau,work,lwork,info)
            !   Form the leading m columns of Q explicitly
            call dorgqr(m,m,m,A,lda,tau,work,lwork,info)

            QQ = A
            RR = matmul(transpose(QQ),AA)
            A  = matmul(QQ,RR) - AA
            err= sum(A)/(m*m*1.0)
            if(err.gt.eps) then
                write(*,*) "Accuracy does not meet the requirements"
                write(*,*) "error:",err,"contrast with eps:",eps
                stop "check point QRsquare_r8 001"
            endif
            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine QRsquare_r8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine PRINT_MATRIX_r4(name,m,A)
    !!=======================================================================76
    !!square matrix A(m,m)
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            character*(*)               name
            integer,intent(in)     ::   m
            real(4),intent(in)     ::   A(m,m)

            !----------------------------------------------------------

            integer             ::  i,  j,  k

            !----------------------------------------------------------

            write(*,*) trim(name)
            do i=1,m
                write(*,"(100(F15.8))") (A(i,j),j=1,m)
            enddo
            write(*,*) "------------------------------------"

            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine PRINT_MATRIX_r4


    subroutine PRINT_MATRIX_r8(name,m,A)
    !!=======================================================================76
    !!square matrix A(m,m)
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            character*(*)               name
            integer,intent(in)     ::   m
            real(8),intent(in)     ::   A(m,m)

            !----------------------------------------------------------

            integer             ::  i,  j,  k

            !----------------------------------------------------------

            write(*,*) trim(name)
            do i=1,m
                write(*,"(100(F15.8))") (A(i,j),j=1,m)
            enddo
            write(*,*) "------------------------------------"

            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine PRINT_MATRIX_r8
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
end module Mathmatic
