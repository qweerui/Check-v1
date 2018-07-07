program main
!!=================================================================================
!!
!!    Tan-Tracker main code
!!
!!    Record of revision:
!!           Date:          Programmer:          Description of change:
!!     ------------       --------------        ------------------------
!!     2018.02.26               Hanrui             NV2.0 original code
!!
!!=================================================================================

    !   common var
    use Common
    use InputOption
    !   Time & Dir assignment
    use TimeCounter,        only : day_back
    use TimeCounter,        only : day_forward
    use TimeCounter,        only : day_counter
    use TimeCounter,        only : SetWinINFO
    !   extract InstCO2
    use extract_bpch,       only : Extract_InstCO2
    use extract_bpch,       only : tmp_CO2_r4
    use extract_bpch,       only : tmp_Prue_r4
    !   reshape
    use reshapeTT,          only : reshape2Dto1D
    use reshapeTT,          only : reshape1Dto2D
    use reshapeTT,          only : reshape3Dto1D
    use reshapeTT,          only : reshape1Dto3D
    !   Obs Obsoperator
    use operator_ocogst,    only : OCOGST_Filename_length
    use operator_ocogst,    only : OCOGST_read
    use operator_ocogst,    only : OCOGST_QualityCtl
    use operator_ocogst,    only : OCOGST_Operator
    !   NLS5 
    use NLSCORE,            only : NLS5_4DVarCORE
    !   statistic
    use Statistic,          only : PercentsLocation_r8
    !   DAInOut  
    use DAInOut,            only : CoeRe_RangeCtl
    use DAInOut,            only : CoeRe_input
    use DAInOut,            only : CoeRe_update
    use DAInOut,            only : CoeRe_errupdate
    use DAInOut,            only : CO2_input
    use DAInOut,            only : CO2_update
    use DAInOut,            only : Flx_extract

    implicit none 
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !                       Time & Dir assignment                     !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        type(WindowsInfo)   ::  NowWin

        !   cycle
        integer             ::  ii, jj, kk, ll
        real                ::  CPU_T1, CPU_T2

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !                         Obsoperator                             !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !   input and output file name
        character(120)      ::  filename_in
        !   reading OCO2
        integer                             ::  OCO2len               !total num day
        type(Satellite_obs_bmc),allocatable ::  daily_OCO2(:)
        !   total obs in current window
        integer                             ::  Obs_length            !total num window
        integer                             ::  Obs_len_lt            !last day num
        integer                             ::  daily_len(Winlength)  !num daily
        type(Satellite_obs_bmc),allocatable ::  Win_obs(:)
        type(Satellite_obs_bmc),allocatable ::  Win_obs_tmp(:)
        !   MPs OPs and Obs after obsoperator
        type(DA_obs_bmc),allocatable        ::  DA_obs(:)

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !                                CO2                              !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !   CO2 DA inputs of NLSCORE
        real(8)             ::  AACO2(nx*ny*nz,NRENS)
        real(8)             ::  X_bCO2(nx*ny*nz)
        real(8)             ::  X_aCO2(nx*ny*nz)
        !   output
        real(8)             ::  X_aErrCO2(nx*ny*nz)
        real(8)             ::  X_a3DCO2(nx,ny,nz)
        real(8)             ::  X_aErr3DCO2(nx,ny,nz)

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !                                Flx                              !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !   FLUX DA inputs of NLSCORE
        !
        !   VAR      NAME                           NTRACER         NLambda     
        !   CO2ff    CO2 fossil fuel emiss                1               1 (cot)
        !   CO2oc    CO2 ocean emissions                  2               3 (sea)
        !   CO2bal   CO2 balanced biosphere               3               2 (lai) (casa)
        !   CO2bb    CO2 biomass burning emiss            4           4 NaN .
        !   CO2bf    CO2 biofuel emission                 5           4 NaN .
        !   CO2nte   CO2 net terr exchange                6               2 (lai)
        !   CO2shp   CO2 ship emissions                   7           4 NaN .
        !   CO2pln   CO2 aircraft emissions               8(3D)       4 NaN .
        !   CO2che   CO2 chemical oxidation               9(3D)       4 NaN .
        !   CO2sur   CO2 chem surf correction            10           4 NaN (oth)
        real(8)             ::  AAFlx(nx*ny*nCoeRe,NRENS)
        real(8)             ::  X_bFlx(nx*ny*nCoeRe)
        real(8)             ::  X_aFlx(nx*ny*nCoeRe)
        !   output
        real(8)             ::  X_aErrFlx(nx*ny*nCoeRe)
        real(8)             ::  X_a3DFlx(nx,ny,nCoeRe)
        real(8)             ::  X_aErr3DFlx(nx,ny,nCoeRe)

        !------------------------------------------------------------------

    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        !                   begin Time & Dir assignment                   !
        !------------------------------------------------------------------
        !   read in nml
        write(*,*) "START Time and Dir assignment"
        call LoadOptions()

        !   set windows info
        NowWin%WCYCLE = NOWCYCLE
        call SetWinINFO(NowWin)
        ! write(*,*) "",NowWin%WCYCLE
        ! write(*,*) "",NowWin%WinStyear
        ! write(*,*) "",NowWin%WinStmonth
        ! write(*,*) "",NowWin%WinStday
        ! write(*,*) "",NowWin%IFstart
        ! write(*,*) "",NowWin%IFend

        !------------------------------------------------------------------
        !                      read in OBS and info                       !
        !------------------------------------------------------------------
        !  OCO2
        write(*,*) "START READING OBS and INFO"
        daily_len  = 0
        Obs_length = 0

        !------------------------------------------------------------------
        !   read oco2 or gosat day by day
        !   20180424 lag Window test RUIHAN
        !   do ii=1, Winlength
        do ii=Foldlen + 1, Winlength

            !   current time
            call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                &ii,WinInyear,WinInmonth,WinInday)
            write(windowinnerdate,'(I8)')  &
                & WinInyear* 10000 + WinInmonth * 100 + WinInday

            !   oco2 file name and length
            ! filename_in = trim(TTDATADIR)
            filename_in = "/data/sciarr/ruihan/OSSEs_TT/OCO2_pseudo/Obs05fix"
            call OCOGST_Filename_length( windowinnerdate, filename_in, OCO2len, "oco")
            !   file exist          length = length
            !   file don't exist    length = -1

            !----------------------------------------------------------
            if( OCO2len.gt.0 ) then

                !--------------------------------------------------
                !   all Obs in one day. allocate daily_OCO2
                allocate(daily_OCO2(OCO2len))

                !   read in obs
                call OCOGST_read( filename_in, OCO2len, daily_OCO2)
                write(*,*) "Latitude range:"
                call PercentsLocation_r8( OCO2len, daily_OCO2(:)%lat*1.0d0, 1.0)

                !   quality control
                filename_in = trim(TTRunDIR)//"/Backrun"
                call OCOGST_QualityCtl( windowinnerdate, filename_in, &
                    &OCO2len, daily_OCO2, "oco")

                !   daily length and windows length
                daily_len(ii) = OCO2len - sum(daily_OCO2(:)%reject)
                Obs_length = sum(daily_len(1:ii))
                Obs_len_lt = Obs_length - daily_len(ii)

                !--------------------------------------------------
                !   SAVE to Win_obs
                !--------------------------------------------------
                !   next part is simple and clever ***
                !   get windows Qctl and thinning obs
                if(Obs_len_lt.ne.0) then
                    allocate(Win_obs_tmp(Obs_len_lt))
                    Win_obs_tmp = Win_obs
                    deallocate(Win_obs)
                endif

                !   total obs in iith day
                allocate(Win_obs(Obs_length))
                kk = 0
                do jj=1,OCO2len
                    if(daily_OCO2(jj)%reject.eq.0) then
                    kk = kk + 1
                    Win_obs( Obs_len_lt + kk ) = daily_OCO2(jj)
                    endif
                enddo

                !   obs in 1:ii-1 th days
                if(Obs_len_lt.ne.0) then
                    Win_obs(1:Obs_len_lt) = Win_obs_tmp(:)
                    deallocate(Win_obs_tmp)
                endif

                !check
                write(*,*) kk,daily_len(ii)
                write(*,*) "Above is info of",ii,"th day in this window."
                write(*,*) ""

                !clean daily_OCO2
                deallocate(daily_OCO2)
                !--------------------------------------------------

            endif
            !----------------------------------------------------------
            write(*,*) "Latitude range:"
            call PercentsLocation_r8( daily_len(ii), &
                & Win_obs(Obs_len_lt+1:Obs_len_lt+daily_len(ii))%lat*1.0d0, 1.0)

        enddo
        !------------------------------------------------------------------
        !   till now we got : Obs_length / daily_len / Win_obs(Obs_length)
        write(*,*) "total obs in DA :", Obs_length, sum(daily_len(:))
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------

    if( Obs_length .gt. 100 ) then

        write(*,*) " Have enough OBS ! "
        !------------------------------------------------------------------
        !                           begin Obsoperator                     !
        !------------------------------------------------------------------

        allocate(DA_obs(Obs_length))

        !------------------------------------------------------------------
        !   do oco2 or gst operator day by day
        call CPU_time(CPU_T1)
        !   20180424 lag Window test RUIHAN
        !   do ii=1, Winlength
        do ii=Foldlen + 1, Winlength

            !   current day info
            if(daily_len(ii).eq.0) cycle
            Obs_len_lt = sum(daily_len(1:ii)) - daily_len(ii)

            !   current time
            call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                &ii,WinInyear,WinInmonth,WinInday)
            write(windowinnerdate,'(I8)')  &
                & WinInyear* 10000 + WinInmonth * 100 + WinInday

            !   obsoperator
            filename_in = trim(TTRunDIR)
            call OCOGST_Operator( windowinnerdate, filename_in, daily_len(ii), &
                & Win_obs(Obs_len_lt+1:Obs_len_lt+daily_len(ii)),  &
                &  DA_obs(Obs_len_lt+1:Obs_len_lt+daily_len(ii)) )

        enddo

        !   time
        call CPU_time(CPU_T2)
        write(*,*) "Time:   ",CPU_T2 - CPU_T1

        !------------------------------------------------------------------

        !------------------------------------------------------------------
        !                   CO2 DA begin the NLS5                         !
        !------------------------------------------------------------------
        if( CO2ENDTimeDA_IOswitch ) then

            !   prepare Xb and AA
            call CO2_input(NowWin,X_bCO2,AACO2)

            !   DA
            write(*,*) "START DA. CO2"
            call NLS5_4DVarCORE( nx, ny, nz, NRENS, &
                    &  Obs_length, DA_obs, AACO2, X_bCO2, X_aCO2, X_aErrCO2)
            write(*,*) "Finish DA. CO2"

            !   output
            call reshape1Dto3D( nx, ny, nz, X_aCO2, X_a3DCO2)
            call reshape1Dto3D( nx, ny, nz, X_aErrCO2, X_aErr3DCO2)
            call CO2_update( NowWin, X_a3DCO2, X_aErr3DCO2)

        endif

        !------------------------------------------------------------------
        !                  FLUX DA begin the NLS5                         !
        !------------------------------------------------------------------
        !   prepare Xb and AA
        call CoeRe_input(NowWin,X_bFlx,AAFlx)

        !   DA
        write(*,*) "START DA. Flux"
        call NLS5_4DVarCORE( nx, ny, nCoeRe, NRENS, &
                &  Obs_length, DA_obs, AAFlx, X_bFlx, X_aFlx, X_aErrFlx)
        write(*,*) "Finish DA. Flux"

        !   time
        call CPU_time(CPU_T1)
        write(*,*) "Time:   ",CPU_T1 - CPU_T2

        !   output
        call reshape3Dto1D( nx, ny, nCoeRe, X_aFlx, X_a3DFlx )
        call reshape3Dto1D( nx, ny, nCoeRe, X_aErrFlx, X_aErr3DFlx )

        !------------------------------------------------------------------
        !             Flux range control, Update and output               !
        !------------------------------------------------------------------
        !   flux range ctl( this win )
        ! call CoeRe_RangeCtl( X_a3DFlx)
        !   flux Update( this win next win )
        call CoeRe_update( NowWin, X_a3DFlx )
        !   flux error output( this win )
        call CoeRe_errupdate( NowWin, X_aErr3DFlx )

        !------------------------------------------------------------------
        !                           cleanup                               !
        !------------------------------------------------------------------
        deallocate(DA_obs)
        deallocate(Win_obs)
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
    else

        write(*,*) " Dont have enough OBS ! Cycle Only "

        !   CO2 Update
        if( CO2ENDTimeDA_IOswitch ) then
            !   prepare Xb and AA
            call CO2_input(NowWin,X_bCO2,AACO2)
            !   output
            X_aCO2 = X_bCO2
            call reshape1Dto3D( nx, ny, nz, X_aCO2, X_a3DCO2)
            X_aErrCO2 = 0.0d0
            call reshape1Dto3D( nx, ny, nz, X_aErrCO2, X_aErr3DCO2)
            call CO2_update( NowWin, X_a3DCO2, X_aErr3DCO2)
        endif

        !   FLUX Update
        !   prepare Xb and AA
        call CoeRe_input(NowWin,X_bFlx,AAFlx)
        !   output
        X_aFlx = X_bFlx
        call reshape3Dto1D( nx, ny, nCoeRe, X_aFlx, X_a3DFlx )
        X_aErrFlx = 0.0d0
        call reshape3Dto1D( nx, ny, nCoeRe, X_aErrFlx, X_aErr3DFlx )
        !   flux Update( this win next win )
        call CoeRe_update( NowWin, X_a3DFlx )
        !   flux error output( this win )
        call CoeRe_errupdate( NowWin, X_aErr3DFlx )

    endif

        !------------------------------------------------------------------
        !                   extract background flux                       !
        !------------------------------------------------------------------
        !   canvert flux into nc
        !   units: g C/m2/day
        call Flx_extract(NowWin, X_a3DFlx, X_aErr3DFlx)
        
        !------------------------------------------------------------------

    !--------------------------------------------------------------------------

!==================================================================================
end program main

