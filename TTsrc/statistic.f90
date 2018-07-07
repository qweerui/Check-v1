MODULE Statistic
!!===============================================================================84
!!
!!  Subroutines of Statistic.
!!
!!    Record of revision:
!!           Date:          Programmer:          Description of change:
!!     ------------       --------------        ------------------------
!!     2014.08.26               Hanrui                   original code
!!     2016.08.04               Hanrui                      integrated
!!     2018.03.28               Hanrui             NV2.0 Final version
!!
!!===============================================================================84
    implicit none
    !--------------------------------------------------------------------------
    PUBLIC          ::  mean_error
                    !    mean error
    PUBLIC          ::  RMS
                    !   Root-mean-square
    PUBLIC          ::  XSigma
                    !   standard Deviation
    PUBLIC          ::  R_T
                    !   correlation coefficient and t-test
    PUBLIC          ::  normalization_Z
                    !   zero-mean normalization_Z
    PUBLIC          ::  normalization_mm
                    !   max_val min_val normalization
    PUBLIC          ::  standardization
                    !   standardization
    PUBLIC          ::  QsortC_r8
                    !   sorts real numbers into ascending numerical order
    PRIVATE         ::  Partition_r8
                    !   QsortC_r8 inner funciton
    PUBLIC          ::  PercentsLocation_r8
                    !   maxval and min val Percents Location of pct
    PUBLIC          ::  CheckDistribution
                    !   Check Distribution

    !--------------------------------------------------------------------------
    CONTAINS


    subroutine mean_error(n,x,base,err)
    !--------------------------------------------------------------------------
    !   this subroutine is used to cacluate the mean error between two
    !   time series 
    !                                                 rui 2014/8/26 
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8),intent(in)          ::  x(n)
        real(8),intent(in)          ::  base(n)
        real(8),intent(out)         ::  err
        !
        real(8)                     ::  temp(n)
        !------------------------------------------------------------------
        temp = 0.0d0
        temp(:) = x(:) - base(:)
        err = sum(temp(:))
        err = err / (n * 1.0d0)
        return
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end  subroutine mean_error


    subroutine RMS(n,x,base,rms_out)
    !--------------------------------------------------------------------------
    !   this subroutine is used to cacluate the root-mean-suqare 
    !   between two time series
    !                                                 rui 2014/8/26 
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8),intent(in)          ::  x(n)
        real(8),intent(in)          ::  base(n)
        real(8),intent(out)         ::  rms_out
        !
        real(8)                     ::  temp_1(n)
        real(8)                     ::  temp(n)
        !------------------------------------------------------------------
        temp_1 = 0.0d0
        temp = 0.0d0
        temp_1(:) = x(:) - base(:)
        !call normalization_mm(n,temp_1,temp)
        temp    = x(:) - base(:)
        temp(:) = temp(:) ** 2
        rms_out = sum(temp(:))
        rms_out = rms_out / (n * 1.0d0)
        rms_out = sqrt(rms_out)
        return
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end  subroutine RMS


    function XSigma(n,x)
    !--------------------------------------------------------------------------
    !   this function is used to cacluate the standard deviation 
    !   of a series
    !       STDEV (sample)       1/(n-1)   √
    !       STDEVP(population)   1/n
    !                                                 rui 2014/8/26 
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8)                     ::  x(n)
        real(8)                     ::  Xsigma
        !
        real(8)                     ::  temp(n),    ave
        !------------------------------------------------------------------

        temp = 0.0d0
        ave = sum(x(:))
        ave = ave / (n * 1.0d0)
        temp(:) = x(:) - ave
        temp(:) = temp(:) ** 2
        Xsigma = sum(temp(:))
        Xsigma = Xsigma / ((n - 1) * 1.0d0)
        Xsigma = sqrt(Xsigma)
        return
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end function XSigma


    subroutine R_T(n,x,y,r_out,t_out)
    !--------------------------------------------------------------------------
    !   this subroutine is used to cacluate the  correlation 
    !   coefficient of two series
    !                                                 rui 2014/8/25 
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)            ::  n
        real(8),intent(in)            ::  x(n)
        real(8),intent(in)            ::  y(n)
        real(8),intent(out)           ::  r_out
        real(8),intent(out)           ::  t_out
        !
        integer             ::  i
        real(8)             ::  ave_x,      ave_y
        real(8)             ::  temp_x(n),  temp_y(n),  multi(n)
        real(8)             ::  numerator,  denominator
        !------------------------------------------------------------------
        !numerator
        ave_x = sum(x(:))
        ave_x = ave_x / (n * 1.0d0)
        ave_y = sum(y(:))
        ave_y = ave_y / (n * 1.0d0)
        temp_x = 0.0d0
        temp_y = 0.0d0
        do i =1,n
        temp_x(i) = x(i) - ave_x
        temp_y(i) = y(i) - ave_y
        multi(i)  = temp_x(i) * temp_y(i)
        enddo
        numerator = sum(multi(:))
        !denominator
        do i =1,n
        temp_x(i) = temp_x(i) ** 2
        temp_y(i) = temp_y(i) ** 2
        enddo
        denominator = sum(temp_x(:))
        denominator = denominator * sum(temp_y(:))
        denominator = sqrt(denominator)
        !R
        r_out = numerator / denominator!*100.0d0
        !T
        t_out = sqrt(1.0d0 - r_out**2)
        t_out = r_out / t_out
        t_out = t_out * sqrt(n*1.0d0 - 2.0d0)
        t_out = abs(t_out)
        return
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine R_T


    subroutine normalization_Z(n,x,x_n)
    !--------------------------------------------------------------------------
    !   zero-mean normalization_Z
    !                                                 rui 2014/12/19
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8),intent(in)          ::  x(n)
        real(8),intent(out)         ::  x_n(n)
        !
        real(8)         ::  temp(n),    ave,      std
        integer         ::  i
        !------------------------------------------------------------------
        ave = sum(x)/(n*1.0d0)
        temp(:) = x(:) - ave
        temp(:) = temp(:)**2
        std = sqrt(sum(temp)/(n*1.0d0))
        do i =1,n
          x_n(i)=(x(i)-ave)/std 
        enddo
        return 
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine normalization_Z


    subroutine normalization_mm(n,x,x_n)
    !--------------------------------------------------------------------------
    !   max_val min_val normalization
    !                                                 rui 2014/12/19
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8),intent(in)          ::  x(n)
        real(8),intent(out)         ::  x_n(n)
        !-------------------------
        real(8)         ::  max_v,      min_v
        integer         ::  i
        !------------------------------------------------------------------
        max_v = maxval(x)
        min_v = minval(x)
        do i =1,n
          x_n(i)=(x(i)-min_v)/(max_v-min_v) 
        enddo
        return 
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine normalization_mm


    subroutine standardization(n,x,x_n)
    !--------------------------------------------------------------------------
    !   standardization
    !                                                 rui 2014/12/19
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        implicit none
        !in & out
        integer,intent(in)          ::  n
        real(8),intent(in)          ::  x(n)
        real(8),intent(out)         ::  x_n(n)
        !
        real(8)         ::  temp(n),    ave,      std
        integer         ::  i
        !------------------------------------------------------------------
        ave = sum(x)/(n*1.0d0)
        x_n(:) = x(:) - ave
        return 
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine standardization

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine QsortC_r8(A)
    !--------------------------------------------------------------------------
    !
    !   download from : http://fcode.cn/code_prof-38-1.html
    !   www.fcode.cn
    !   Recursive Fortran 95 quicksort routine
    !   sorts real numbers into ascending numerical order
    !   Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
    !   Based on algorithm from Cormen et al., Introduction to Algorithms,
    !   1997 printing
    !   Made F conformant by Walt Brainerd
    !   
    !--------------------------------------------------------------------------

        !-----------------------------------------------------------------
        real(8), intent(in out), dimension(:) :: A
        integer :: iq

        !-----------------------------------------------------------------
        if(size(A) > 1) then
            call Partition_r8(A, iq)
            call QsortC_r8(A(:iq-1))
            call QsortC_r8(A(iq:))
        endif

        !-----------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine QsortC_r8

    subroutine Partition_r8(A, marker)
    !--------------------------------------------------------------------------
    !   download from : http://fcode.cn/code_prof-38-1.html
    !   www.fcode.cn
    !--------------------------------------------------------------------------

        !-----------------------------------------------------------------
        real(8), intent(in out), dimension(:) :: A
        integer, intent(out) :: marker
        integer     :: i, j
        real(8)     :: temp
        real(8)     :: x      ! pivot point

        !-----------------------------------------------------------------
        x = A(1)
        i = 0
        j = size(A) + 1
        do
        j = j-1
            do
            if (A(j) <= x) exit
            j = j-1
            end do
        i = i+1
            do
            if (A(i) >= x) exit
            i = i+1
            end do
        if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
        elseif (i == j) then
            marker = i+1
            return
        else
            marker = i
        return
        endif
        end do
        !-----------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine Partition_r8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine PercentsLocation_r8(dim,AA,pct)
    !--------------------------------------------------------------------------
    !   maxval and min val Percents Location of pct
    !   pct 0.0d0 ~ 1.0d0
    !--------------------------------------------------------------------------
    
        !-----------------------------------------------------------------

        implicit none
        integer         ::  dim
        real(8)         ::  AA(dim)
        real            ::  pct
        real(8)         ::  AAmin,AAmax
        !-----------------------------------------------------------------
        real(8)         ::  AArank(dim)
        integer         ::  i,j,k
        integer         ::  minum,maxum
        !-----------------------------------------------------------------
        AArank(:) = AA(:)
        call QsortC_r8(AArank)
        minum = int(dim*(1-pct))
        if(minum.eq.0) minum = 1
        AAmin=AArank(minum)
        maxum = int(dim*pct)
        AAmax=AArank(maxum)
        write(*,*) "-   -   -   -   -   -   -   -   -   -   "
        write(*,*) "Percentage :   ",int(pct*100),"%"
        write(*,*) "Max Value  :   ",AAmax
        write(*,*) "Min Value  :   ",AAmin
        write(*,*) ""
        return
        !-----------------------------------------------------------------

    !--------------------------------------------------------------------------
    end subroutine PercentsLocation_r8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CheckDistribution(num,array,NN,ObsReject)
    !!=======================================================================76
    !!
    !!  check the distribution of error between OBS and SIM
    !!  flitered the obs lied out of 3Sigma(2Sigma)
    !!
    !!=======================================================================76
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input info
            integer,intent(in)  ::  num
            real(8),intent(in)  ::  array(num)
            integer,intent(in)  ::  NN
            integer,intent(out) ::  ObsReject(num)
            !   array info
            real(8)             ::  average
            real(8)             ::  Sigma
            !   box info
            real(8)             ::  startpoint
            real(8)             ::  endpoint
            real(8)             ::  boxwidth
            integer             ::  boxnum
            integer,allocatable ::  box(:)
            real(8),allocatable ::  boxcenter(:)
            integer             ::  serialnumber 
            !   parameters
            integer             ::  i,  j,  k
            integer             ::  ii, jj, kk
            integer             ::  unit1

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   basic info of X
            unit1 = 111
            open(unit1,file="SigmaCheck")
            write(unit1,*) "Check Distribution : "
            write(unit1,*) "-----------------------------------------"
            write(unit1,*) "maxval value of X:  = ",maxval(array)
            write(unit1,*) "minval value of X:  = ",minval(array)
            average = sum(array)/(1.0d0*num)
            write(unit1,*) "average value of X: = ",average
            Sigma = XSigma(num,array)
            write(unit1,*) "Std value of X:     = ",Sigma

            !   2~3Sigma start and end
            startpoint = floor(average - NN*Sigma)
            endpoint = ceiling(average + NN*Sigma)
            write(unit1,*) "Start and End point:"
            write(unit1,"(3F10.2)") startpoint,endpoint,endpoint-startpoint
            write(unit1,*) "-----------------------------------------"

            !   box width and box num
            if(endpoint-startpoint .lt.  0.0d0) then
                write(*,*) "endpoint .lt. startpoint sub CheckDistribution 01"
                stop
            endif
            if(endpoint-startpoint .le.  1.0d0) boxwidth = 0.025d0
            if(endpoint-startpoint .gt.  1.0d0) boxwidth = 0.05d0
            if(endpoint-startpoint .gt.  2.0d0) boxwidth = 0.1d0
            if(endpoint-startpoint .gt.  5.0d0) boxwidth = 0.2d0
            if(endpoint-startpoint .gt. 10.0d0) boxwidth = 0.5d0
            if(endpoint-startpoint .gt. 20.0d0) boxwidth = 1.0d0
            if(endpoint-startpoint .gt. 40.0d0) boxwidth = 2.0d0
            boxwidth = boxwidth/2.0d0
            boxnum = int((endpoint - startpoint)/boxwidth) + 2

            !   allocate
            allocate(box(boxnum))
            allocate(boxcenter(boxnum))

            !   boxcenter
            do jj=1,boxnum
                boxcenter(jj) = startpoint + (jj-1)*boxwidth - boxwidth/2.0d0 
            enddo

            !   filter
            box = 0
            ObsReject = 0
            do ii=1,num
                if(array(ii).lt.startpoint) then
                box(1) = box(1) + 1
                ObsReject(ii) = 1
                cycle
                endif
                if(array(ii).gt.endpoint) then
                box(boxnum) = box(boxnum) + 1
                ObsReject(ii) = 1
                cycle
                endif
                serialnumber = floor((array(ii) - startpoint)/boxwidth) + 2
                box(serialnumber) = box(serialnumber) + 1
            enddo

            !   distribution
            write(unit1,*) "box num"
            do jj=1,boxnum
                write(unit1,"(I8,F10.2,I8)") jj,boxcenter(jj),box(jj)
            enddo
            write(unit1,*) "-----------------------------------------"

            !   deallocate
            deallocate(box)

            !----------------------------------------------------------

        !------------------------------------------------------------------
        
    !!=======================================================================76
    end subroutine CheckDistribution

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
end MODULE Statistic
