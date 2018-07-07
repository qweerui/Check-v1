Module Random
!!===============================================================================84
!!
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2016.06.22                fcode                   original code:
!!                                                                     ran
!!                                                                  normal
!!
!!        2016.06.22               Hanrui         Random_Orthogonal_Matrix
!!
!!        2016.08.02               Hanrui                      ImpSampling
!!        2018.03.28               Hanrui              NV2.0 Final version
!!                                                  
!!===============================================================================84

    implicit None
    !--------------------------------------------------------------------------
    PUBLIC          ::  ran               
                    !   return a uniform random number between 0-1 
    PUBLIC          ::  normal   
                    !   return a normal distribution
    PUBLIC          ::  Random_Orthogonal_Matrix
                    !   generate random orthgnonal matrix ROMA(beta*N,beta*N)
    PUBLIC          ::  ImpSampling
                    !   Improved sampling strategy Data Assimilation (geir Evensen)


    !--------------------------------------------------------------------------
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function ran()
    !!=======================================================================76
    !! returns random number between 0 - 1 
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------
            integer , save  :: flag_ranModsave = 0
            real(8)         :: ran

            if(flag_ranModsave.eq.0) then
              call random_seed()
              flag_ranModsave = 1
            endif
            call random_number(ran)
        !------------------------------------------------------------------

    !!=======================================================================76
    end function ran
   
    function normal(mean,sigma)
    !!=======================================================================76
    !!  return random normal distributed number with 
    !!  mean value: mean and standard deviation: sigma
    !!=======================================================================76
    !!The above codes are made in Fortran 90 language
    !!if you have any question, you may write to  sealin2008@hotmail.com
    !!
    !!parts of code to test function ran() and subroutine normal(mean,sigma):
    !!mean  = 0
    !!sigma = 1
    !!x~N(0,1)    standarded normal distribution
    !!3sigma principle :
    !!counting numbers in series x distributed in (-3,3).
    !!
    !
    !real*8                  ::  randAA(1000000)
    !integer                 ::  rank(60)
    !
    !rank = 0
    !do i=1,1000000
    !    randAA(i) = normal(0.0d0,1.0d0)
    !    if(floor(randAA(i)*10)+31.ge.1.and.floor(randAA(i)*10)+31.le.60) then
    !        rank(floor(randAA(i)*10)+31) = rank(floor(randAA(i)*10)+31) + 1
    !    endif
    !enddo
    !
    !open(1,file="1.txt")
    !do j=1,60
    !write(1,*) rank(j)
    !enddo
    !close(1)

        implicit none
        !------------------------------------------------------------------

            integer :: flag_ranModsave
            real(8), parameter :: pi = 3.141592653589793239 
            real(8) :: u1, u2, y1, y2, normal, mean, sigma
            save flag_ranModsave
            data flag_ranModsave /0/

        !------------------------------------------------------------------

            u1 = ran()
            u2 = ran()
            if (flag_ranModsave.eq.0) then
                y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
                normal = mean + sigma*y1
                flag_ranModsave = 1
            else
                y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
                normal = mean + sigma*y2
                flag_ranModsave = 0
            endif 

        !------------------------------------------------------------------

    !!=======================================================================76
    end function normal

    function ran_r4()
    !!=======================================================================76
    !! returns random number between 0 - 1 
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------
            integer , save  :: flag_ranModsave = 0
            real            :: ran_r4

            if(flag_ranModsave.eq.0) then
              call random_seed()
              flag_ranModsave = 1
            endif
            call random_number(ran_r4)
        !------------------------------------------------------------------

    !!=======================================================================76
    end function ran_r4
   
    function normal_r4(mean,sigma)
    !!=======================================================================76
    !!  return random normal distributed number with 
    !!  mean value: mean and standard deviation: sigma
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            integer :: flag_ranModsave
            real   , parameter :: pi = 3.141592653589793239 
            real    :: u1, u2, y1, y2, normal_r4, mean, sigma
            save flag_ranModsave
            data flag_ranModsave /0/

        !------------------------------------------------------------------

            u1 = ran_r4()
            u2 = ran_r4()
            if (flag_ranModsave.eq.0) then
                y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
                normal_r4 = mean + sigma*y1
                flag_ranModsave = 1
            else
                y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
                normal_r4 = mean + sigma*y2
                flag_ranModsave = 0
            endif 

        !------------------------------------------------------------------

    !!=======================================================================76
    end function normal_r4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Random_Orthogonal_Matrix(beta,N,ROMA)
    !!=======================================================================76
    !!generate random orthgnonal matrix ROMA(beta*N,beta*N)
    !!
    !!  ROMA(:,i)  .O.  ROMA(:,i)   checked
    !!  1. write(*,*) matmul(transpose(Q(:,1:2)),Q(:,1:2))
    !!  2. leading column of ROMA is the orthogonal vector(lapack manual).
    !!
    !!=======================================================================76

        use Mathmatic,          only :  QRsquare_r8
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            integer,intent(in)      ::  beta
            integer,intent(in)      ::  N
            real(8),intent(inout)   ::  ROMA(beta*N,beta*N)

            !----------------------------------------------------------
            !   parameters
            integer                 ::  i,  j,  k

            !----------------------------------------------------------
            !   identity matrix I  
            real(8)                 ::  II(beta*N,beta*N)
            !   Q,R
            real(8)                 ::  QQ(beta*N,beta*N)
            real(8)                 ::  RR(beta*N,beta*N)

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------

            do i = 1,beta*N
            do j = 1,beta*N
                ROMA(i,j) = normal(0.0d0,1.0d0)
            enddo
            enddo


            call QRsquare_r8(beta*N,ROMA,QQ,RR)

            II = 0.0
            do i=1,beta*N
                if(RR(i,i).gt.0.0) then
                    II(i,i) =  1.0
                else if(RR(i,i).lt.0.0) then
                    II(i,i) = -1.0
                else
                    stop "failed in QR factorization 002"
                endif
            enddo

            !!  checked!
            !!  call PRINT_MATRIX_r8("Q",m,QQ)
            !!  call PRINT_MATRIX_r8("R",m,RR)
            !!  tmp = matmul(transpose(QQ),QQ)
            !!  call PRINT_MATRIX_r8("QTQ",m,tmp)
            !!  tmp = matmul(QQ,RR) - ROMA
            !!  call PRINT_MATRIX_r8("QR-A",m,tmp)

            QQ = matmul(QQ,II)
            ROMA = QQ

            !----------------------------------------------------------

        !------------------------------------------------------------------
        
    !!=======================================================================76
    end subroutine Random_Orthogonal_Matrix

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ImpSampling(ndim,alpha,beta,nrens,AA,AAnew)
    !!=======================================================================76
    !!
    !!  Improved sampling strategy page 163 Data Assimilation (geir Evensen)
    !!
    !!=======================================================================76

        use Mathmatic,          only : SVD_BasedOnEOF_r8
        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   variable parameters
            integer,intent(in)      ::  ndim
            integer,intent(in)      ::  alpha
            integer,intent(in)      ::  beta
            integer,intent(in)      ::  nrens

            !   large ensemble singular value decomposition:
            real(8),intent(in)      ::  AA(ndim,alpha*nrens)
            real(8)                 ::  AABak(ndim,alpha*nrens)
            real(8)                 ::  UU(ndim,alpha*nrens)
            real(8)                 ::  SigMa(alpha*nrens,alpha*nrens)
            real(8)                 ::  VVT(alpha*nrens,alpha*nrens)

            !   new ensemble sampled:
            real(8),intent(inout)   ::  AAnew(ndim,nrens)
            real(8)                 ::  UUR(ndim,beta*nrens)
            real(8)                 ::  SigMaR(beta*nrens,beta*nrens)
            real(8)                 ::  VVTR(beta*nrens,beta*nrens)
            real(8)                 ::  VVTR_Extract(beta*nrens,nrens)

            !   parameter
            integer                 ::  i,  j,  k

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   singular value decomposition
            call SVD_BasedOnEOF_r8(ndim,alpha*nrens,AA,UU,SigMa,VVT)

            !   extract UU to UUR
            !   extract SigMa to SigMaR
            SigMaR = 0.0
            do i=1,beta*nrens
                UUR(:,i) = UU(:,i)
                SigMaR(i,i) = sqrt((beta*1.0)/(alpha*1.0)) * SigMa(i,i)
            enddo

            !   built random orthogonal matrix VVTR
            call Random_Orthogonal_Matrix(beta,nrens,VVTR)
            do i=1,nrens
                VVTR_Extract(:,i) = VVTR(:,i)
            enddo

            !   reform new ensemble AA
            AAnew = matmul(matmul(UUR,SigMaR),VVTR_Extract)

            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine ImpSampling

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
End Module Random
