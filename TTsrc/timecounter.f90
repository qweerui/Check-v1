module TimeCounter
!!===============================================================================84
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2014.10.08               Hanrui           original code

      !counting strategy :
                                     !=======!
                                     !   T   !
      !---------------------------------------------------------------------!
      !             1        2  3  4  5  6  7  8  9        10               !
      !          2010.1.1      .  .  .  .  .  .        2010.1.10            !
      !          total namber = 10                                          !
      !---------------------------------------------------------------------!

                                     !=======!
                                     !   F   !
      !---------------------------------------------------------------------!
      !                                                                     !
      !             0        1  2  3  4  5  6  7  8        9                !
      !          2010.1.1      .  .  .  .  .  .        2010.1.10            !
      !          total namber = 9                                           !
      !                                                                     !
      !             1        2  3  4  5  6  7  8  9        0                !
      !          2010.1.1      .  .  .  .  .  .        2010.1.10            !
      !          total namber = 9                                           !
      !                                                                     !
      !             0       1  2  3  4  5  6  7  8         0                !
      !          2010.1.1      .  .  .  .  .  .        2010.1.10            !
      !          total namber = 8                                           !
      !                                                                     !
      !=====================================================================!

!!    Record of revision:
!!           Date:          Programmer:          Description of change:
!!     ------------       --------------        ------------------------
!!     2017.10.29               Hanrui             NV2.0 original code
!!     2017.10.29               Hanrui                      SetWinINFO
!!     2018.03.28               Hanrui             NV2.0 Final version
!!
!!===============================================================================84

    implicit none
    !--------------------------------------------------------------------------
    PUBLIC          ::  day_forward
                    !   To give the time_num th date from the start date
    PUBLIC          ::  day_counter
                    !   To count days between start date to end date
    PUBLIC          ::  day_back
                    !   To give the start date from the time_num th date
    PRIVATE         ::  leapyear
                    !   To find out whether a year is leap or not
    PUBLIC          ::  SetWinINFO
                    !   Give full window time info from cyclenum

    !--------------------------------------------------------------------------
    contains

      subroutine day_forward(styear,stmon,stday,time_num,year,month,day)
      !=====================================================================
      !this subroutine is used to transform time series number into date
      !rui  2014/8/13
      !                          same effect with :AssTime
      !                                            Tian Xiangjun,2013/02/18
      !=====================================================================

        use Common,only  :  Monthdays
        implicit none

        !---------------------------------------------------------------
        !year,month,day        :  current
        !time_num              :  current numerical order in time series
        !styear,stmon,stday  :  start 
        !---------------------------------------------------------------
  
        integer, intent(in)     :: styear,  stmon,  stday,  time_num
        integer, intent(inout)  :: year,    month,  day
        logical                 :: leap
        integer                 :: Mdays(12)

        !---------------------------------------------------------------
        !---original data---
        day   = stday + time_num - 1
        month = stmon
        year  = styear
        Mdays = Monthdays

          call leapyear(year,leap)
          if(leap) then
            Mdays(2) = 29
          else
            Mdays(2) = 28
          endif
        !---judgement---
        do while(day .gt. Mdays(month))
          if(month .eq. 12) then
            year = year + 1
              call leapyear(year,leap)
              if(leap) then
                Mdays(2) = 29
              else
                Mdays(2) = 28
              endif
            month = 1
            day = day - Mdays(12)
          else
            day = day - Mdays(month)
            month = month + 1
          endif
        enddo
        if(day == 0)  day = Mdays(month)
        !---------------------------------------------------------------
        return
      !=====================================================================
      end subroutine day_forward

      subroutine day_counter(styear,stmon,stday,edyear,edmon,edday,i)
      !=====================================================================
      !count days between start date to end date
      !rui  2014/8/15
      !=====================================================================

        use Common,only  :  Monthdays
        implicit none 

        !---------------------------------------------------------------

        integer, intent(in)       ::  edyear,  edmon,  edday
        integer, intent(in)       ::  styear,  stmon,  stday
        integer, intent(inout)    ::  i
        integer                   ::  k
        logical                   ::  leap
        integer                   ::  Mdays(12)


        !---------------------------------------------------------------
        !days in entire year
        i = 0
          do k = styear,edyear
            call leapyear(k,leap)
              if ( leap ) then 
                i = i + 366
              else
                i = i + 365
              endif
          enddo
        !days in styear and edyear
        Mdays = Monthdays
          call leapyear(styear,leap)
          if(leap)  Mdays(2) = 29
          do k = 1 , stmon
            i = i - Mdays(k)
          enddo
          i = i + Mdays(stmon) - stday

          !give back the original data of Mdays(2)
          Mdays(2) = 28

          call leapyear(edyear,leap)
          if(leap)  Mdays(2) = 29
          do k = edmon , 12
            i = i - Mdays(k) 
          enddo
          i = i + edday + 1

        !---------------------------------------------------------------
        return
      !=====================================================================
      end subroutine day_counter

      subroutine day_back(edyear,edmon,edday,time_num,year,month,day)
      !=====================================================================
      !to find the start date with the date number and end date
      !rui  2014/10/8
      !=====================================================================

        use Common,only  :  Monthdays
        implicit none

        !---------------------------------------------------------------
        !edyear,demon,edday    :  current
        !time_num              :  current numerical order in time series
        !styear,stmon,stday  :  start 
        !---------------------------------------------------------------

        integer, intent(in)     :: edyear,  edmon,  edday,  time_num
        integer, intent(inout)  :: year,    month,    day
        logical                 :: leap
        integer                 :: Mdays(12)

        !---------------------------------------------------------------
        !---original data---
        day   = edday - time_num + 1
        month = edmon
        year  = edyear
        Mdays = Monthdays

          call leapyear(year,leap)
          if(leap) then
            Mdays(2) = 29
          else
            Mdays(2) = 28
          endif
        !---judgement---
        do while(day .le. 0)
          if(month == 1) then
            year = year - 1
              call leapyear(year,leap)
              if(leap) then
                Mdays(2) = 29
              else
                Mdays(2) = 28
              endif
            month = 12
            day = day + Mdays(12)
          else
            day = day + Mdays(month-1)
            month = month - 1
          endif
        enddo

        !---------------------------------------------------------------
        return
      !=====================================================================
      end subroutine day_back

      subroutine leapyear(year,leap)
      !=====================================================================
      !this subroutine is used to find out leap year
      !rui  2014/8/14
      !=====================================================================
        implicit none
        integer, intent(in)     :: year
        logical, intent(inout)  :: leap
        !----------------------------------------------------------------
        !---leap year or not---
        if (mod(year,100)==0) then
          if (mod(year,400)==0) then
          leap = .true.
          endif
        elseif(mod(year,4)==0) then
          leap = .true.
        else 
          leap = .false.
        endif
        return
      end subroutine leapyear
      !=====================================================================

      subroutine SetWinINFO(CycleWin)
      !!====================================================================
      !!  cyclenum -> CycleWin
      !!  rui  2017/10/29
      !!====================================================================
          use Common,         only    :   Monthdays
          use InputOption,    only    :   WindowsInfo
          use InputOption,    only    :   Winlength
          use InputOption,    only    :   Fold_IOswitch,  Foldlen
          use InputOption,    only    :   STARTYEAR,      STARTMONTH
          use InputOption,    only    :   ENDYEAR,        ENDMONTH
          implicit none
          !-------------------------------------------------------------
          type(WindowsInfo),intent(inout) ::  CycleWin

          !-------------------------------------------------------------
          integer             ::  passlength  ! days passed
          integer             ::  Mdays(12)   
          integer             ::  Total_len
          logical             ::  leap
          integer             ::  cyclenum

          !-------------------------------------------------------------
          cyclenum = CycleWin%WCYCLE

          !-------------------------------------------------------------
          ! caculate passed day length
          passlength = 1 + Winlength * (cyclenum - 1)
          if(Fold_IOswitch) then
            passlength = 1 + (Winlength - Foldlen) * (cyclenum - 1)
          endif

          !-------------------------------------------------------------
          ! set start year month day
          call day_forward(STARTYEAR,STARTMONTH,01,passlength,&
              &CycleWin%WinStyear,CycleWin%WinStmonth,CycleWin%WinStday)
          
          CycleWin%WCYCLE = cyclenum
          
          !-------------------------------------------------------------
          ! set IFstart
          CycleWin%IFstart = .false.
          CycleWin%IFend   = .false.
          if(cyclenum.eq.1) CycleWin%IFstart = .true.

          !-------------------------------------------------------------
          ! set IFend
          Mdays = Monthdays
          call leapyear(ENDYEAR,leap)
          if(leap) Mdays(2) = 29
          call day_counter(STARTYEAR,STARTMONTH,01,&
              &ENDYEAR,ENDMONTH,Mdays(ENDMONTH),Total_len)
          if(Total_len.le.passlength) then
            CycleWin%IFend = .true.
          endif

          !-------------------------------------------------------------
          return

      !!====================================================================
      end subroutine SetWinINFO

    !--------------------------------------------------------------------------
    
!!===============================================================================84
end module TimeCounter
