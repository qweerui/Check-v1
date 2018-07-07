MODULE DAInOut
!!=================================================================================
!!
!!   flux Xb AA prepare, flux range ctl, flux update and flux output
!!   CO2 Xb AA prepare and CO2 output
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2017.04.02                   Hanrui                 NV2.0 original code
!!
!!=================================================================================

    !   common var
    use Common
    use InputOption
    !   Flux range ctl
    use Statistic,          only : XSigma
    !   read in CoeRe and CO2
    use readwrite_nc,       only : Read_CoeRe_nc
    use readwrite_nc,       only : tmp_CoeRe_r4
    use readwrite_nc,       only : Write_CoeRe_nc
    use readwrite_nc,       only : Read_CO2_nc
    use readwrite_nc,       only : tmp_CO2nc_r4
    use readwrite_nc,       only : Write_CO2_nc
    !   extract Flux
    use readwrite_nc,       only : tmp_Flux_r4
    use readwrite_nc,       only : Write_Flux_nc
    use extract_bpch,       only : Extract_Flux
    use extract_bpch,       only : tmp_Fluxwinlen_r4
    !   Time & Dir assignment
    use TimeCounter,        only : day_forward
    use TimeCounter,        only : SetWinINFO
    !   reshape
    use reshapeTT,          only : reshape2Dto1D
    use reshapeTT,          only : reshape1Dto2D
    use reshapeTT,          only : reshape3Dto1D
    use reshapeTT,          only : reshape1Dto3D

    IMPLICIT NONE

    !--------------------------------------------------------------------------
    PUBLIC                  ::  CoeRe_RangeCtl
                                ! Flux range ctl
    PUBLIC                  ::  CoeRe_update
                                ! Flux update this win next win
    PUBLIC                  ::  CoeRe_errupdate
                                ! Flux error out put
    PUBLIC                  ::  CoeRe_input
                                ! prepare Xb and AA for FluxDA

    !--------------------------------------------------------------------------
    PUBLIC                  ::  CO2_input
                                ! prepare Xb and AA for CO2DA
    PUBLIC                  ::  CO2_update
                                ! CO2 update this win next win

    !--------------------------------------------------------------------------
    PUBLIC                  ::  Flx_extract
                                ! convert Flux from bpch to nc

    !--------------------------------------------------------------------------
    character(120),PRIVATE  ::  filename_in

    !--------------------------------------------------------------------------
    Contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CoeRe_RangeCtl(X_a3DFlx)
    !!=========================================================================
    !!
    !!  flux range control
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            real(8),intent(inout)   ::  X_a3DFlx(nx,ny,nCoeRe)
            !   xa range ctl
            real(8)                 ::  X_aFlxlev(nx*ny)
            real(8)                 ::  average,    SigMa
            !   loop
            integer                 ::  ii, jj, kk
            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------

            do ii=1,nCoeRe

                !   nCoeRe: 1 2 3 6
                if(ii.eq.1 .or. ii.eq.2 .or. ii.eq.3 .or. ii.eq.6 ) then

                call reshape2Dto1D( nx, ny, X_a3DFlx(:,:,ii), X_aFlxlev )
                write(*,*) "Maxval Flux Xa:",maxval(X_aFlxlev)
                write(*,*) "Minval Flux Xa:",minval(X_aFlxlev)

                ! average = sum(X_aFlxlev)/( nx*ny*1.0d0)
                ! SigMa = XSigma( nx*ny, X_aFlxlev)
                average = 1.0d0
                if(ii.eq.1) then
                    SigMa = SamplescaleCot
                elseif(ii.eq.2) then
                    SigMa = SamplescaleSea
                elseif(ii.eq.3 .or. ii.eq.6) then
                    SigMa = SamplescaleLai
                endif
                write(*,*) "Average Flux Xa:",average
                write(*,*) "Standarded deviation Flux Xa:",SigMa
                write(*,*) "Range Xa:",average - 3*SigMa, average + 3*SigMa

                do jj=1,nx*ny
                if( X_aFlxlev(jj) .gt. average + 3*SigMa ) then
                    X_aFlxlev(jj) = average + 3*SigMa
                elseif( X_aFlxlev(jj) .lt. average - 3*SigMa ) then
                    X_aFlxlev(jj) = average + 3*SigMa
                endif
                enddo
                write(*,*) "Maxval Flux Xa (range ctl):",maxval(X_aFlxlev)
                write(*,*) "Minval Flux Xa (range ctl):",minval(X_aFlxlev)
                call reshape1Dto2D( nx, ny, X_aFlxlev, X_a3DFlx(:,:,ii) )

                endif

            enddo
            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CoeRe_RangeCtl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CoeRe_update(NowWin,X_a3DFlx)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(in)              ::  X_a3DFlx(nx,ny,nCoeRe)

            !   next win and pre win time
            type(WindowsInfo)       ::  PreWin
            type(WindowsInfo)       ::  NxtWin
            !   win Update
            real(8)                 ::  NxtWin_3DFlx(nx,ny,nCoeRe)
            real(8)                 ::  PreWin_3DFlx(nx,ny,nCoeRe)
            !   loop
            integer                 ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   out put (this win)
            tmp_CoeRe_r4(:,:,:,1) = real(X_a3DFlx(:,:,:))
            !   cycle
            do ii=1,Winlength
                !   current time
                call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                    &ii,WinInyear,WinInmonth,WinInday)
                write(windowinnerdate,'(I8)')  &
                    & WinInyear* 10000 + WinInmonth * 100 + WinInday
                filename_in = trim(TTRunDIR)//"/Backrun/CoeRe/CoeRe"//windowinnerdate//".nc"
                call Write_CoeRe_nc( filename_in, windowinnerdate, tmp_CoeRe_r4 )
            enddo

            !----------------------------------------------------------
            !   out put (next win)
            !   prevous win updated flux(Xa_pres)
            if(NowWin%IFstart) then
                tmp_CoeRe_r4 = 1.0
            else
                PreWin%WCYCLE = NowWin%WCYCLE - 1
                call SetWinINFO(PreWin)
                write(windowinnerdate,'(I8)')  &
                    & PreWin%WinStyear* 10000 + PreWin%WinStmonth * 100 + PreWin%WinStday
                filename_in = trim(TTRunDIR)//"/Backrun/CoeRe/CoeRe"//windowinnerdate//".nc"
                call Read_CoeRe_nc( filename_in, tmp_CoeRe_r4)
            endif
            PreWin_3DFlx(:,:,:) = tmp_CoeRe_r4(:,:,:,1) * 1.0d0
            NxtWin_3DFlx(:,:,:) = PreWin_3DFlx(:,:,:) + X_a3DFlx(:,:,:) + 1.0d0
            NxtWin_3DFlx(:,:,:) = NxtWin_3DFlx(:,:,:)/3.0d0

            !   next win updated flux(Xb_nxt)
            NxtWin%WCYCLE = NowWin%WCYCLE  + 1
            call SetWinINFO(NxtWin)
            tmp_CoeRe_r4(:,:,:,1) = real(NxtWin_3DFlx(:,:,:))
            !   cycle
            do ii=1,Winlength
                !   current time
                call day_forward(NxtWin%WinStyear,NxtWin%WinStmonth,NxtWin%WinStday,&
                    &ii,WinInyear,WinInmonth,WinInday)
                write(windowinnerdate,'(I8)')  &
                    & WinInyear* 10000 + WinInmonth * 100 + WinInday
                filename_in = trim(TTRunDIR)//"/Backrun/CoeRe/CoeRe"//windowinnerdate//".nc"
                call Write_CoeRe_nc( filename_in, windowinnerdate, tmp_CoeRe_r4 )
            enddo

            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CoeRe_update

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CoeRe_errupdate(NowWin,X_aErr3DFlx)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(in)              ::  X_aErr3DFlx(nx,ny,nCoeRe)

            !   loop
            integer                         ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   out put err (this win)
            tmp_CoeRe_r4(:,:,:,1) = real(X_aErr3DFlx(:,:,:))
            !   cycle
            do ii=1,Winlength
                !   current time
                call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                    &ii,WinInyear,WinInmonth,WinInday)
                write(windowinnerdate,'(I8)')  &
                    & WinInyear* 10000 + WinInmonth * 100 + WinInday
                filename_in = trim(TTRunDIR)//"/Backrun/CoeRe/CoeReErr"//windowinnerdate//".nc"
                call Write_CoeRe_nc( filename_in, windowinnerdate, tmp_CoeRe_r4 )
            enddo

            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CoeRe_errupdate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CoeRe_input(NowWin,X_bFlx,AAFlx)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(out)             ::  X_bFlx(nx*ny*nCoeRe)
            real(8),intent(out)             ::  AAFlx(nx*ny*nCoeRe,NRENS)

            !   loop
            integer                         ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   analysis time
            write(windowinnerdate,'(I8)') NowWin%WinStyear * 10000 &
            &+ NowWin%WinStmonth * 100 + NowWin%WinStday

            !   Xb
            filename_in = trim(TTRunDIR)//"/Backrun/CoeRe/CoeRe"//windowinnerdate//".nc"
            call Read_CoeRe_nc( filename_in, tmp_CoeRe_r4)
            call reshape3Dto1D(nx,ny,nCoeRe,tmp_CoeRe_r4(:,:,:,1)*1.0d0,X_bFlx)

            !   samples AA
            do jj =1,NRENS
            write(windowinnerdate,'(I8)') NowWin%WinStyear * 10000 &
            &+ NowWin%WinStmonth * 100 + NowWin%WinStday
            filename_in = trim(TTRunDIR)//"/"//NumCha(jj)//"/CoeRe/CoeRe"//windowinnerdate//".nc"
            call Read_CoeRe_nc( filename_in, tmp_CoeRe_r4)
            call reshape3Dto1D(nx,ny,nCoeRe,tmp_CoeRe_r4(:,:,:,1)*1.0d0,AAFlx(:,jj))
            enddo

            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CoeRe_input

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CO2_input(NowWin,X_bCO2,AACO2)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(out)             ::  X_bCO2(nx*ny*nz)
            real(8),intent(out)             ::  AACO2(nx*ny*nz,NRENS)

            !   loop
            integer                         ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   analysis time
            call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                    &Winlength+1,WinInyear,WinInmonth,WinInday)
            write(windowinnerdate,'(I8)') WinInyear* 10000 + WinInmonth * 100 + WinInday

            !   Xb
            filename_in = trim(TTRunDIR)//"/Backrun/GEOSChem_restart."&
                    &//windowinnerdate//"0000.nc"
            call Read_CO2_nc( filename_in, tmp_CO2nc_r4)
            call reshape3Dto1D(nx,ny,nz,tmp_CO2nc_r4(:,:,:,1)*ppm*1.0d0,X_bCO2)

            !   samples AA
            do jj =1,NRENS
            filename_in = trim(TTRunDIR)//"/"//NumCha(jj)//"/GEOSChem_restart."&
                    &//windowinnerdate//"0000.nc"
            call Read_CO2_nc( filename_in, tmp_CO2nc_r4)
            call reshape3Dto1D(nx,ny,nz,tmp_CO2nc_r4(:,:,:,1)*ppm*1.0d0,AACO2(:,jj))
            enddo

            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CO2_input


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CO2_update(NowWin,X_a3DCO2,X_aErr3DCO2)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(in)              ::  X_a3DCO2(nx,ny,nz)
            real(8),intent(in)              ::  X_aErr3DCO2(nx,ny,nz)

            !   loop
            integer                         ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   analysis time
            call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                    &Winlength+1,WinInyear,WinInmonth,WinInday)
            write(windowinnerdate,'(I8)') WinInyear* 10000 + WinInmonth * 100 + WinInday

            !   output co2
            filename_in = trim(TTRunDIR)//"/Backrun/GEOSChem_restart.da."//windowinnerdate//"0000.nc"
            tmp_CO2nc_r4(:,:,:,1) = real(X_a3DCO2(:,:,:)/ppm )
            call Write_CO2_nc( filename_in, windowinnerdate, tmp_CO2nc_r4 )

            !   output co2 err
            filename_in = trim(TTRunDIR)//"/Backrun/GEOSChem_restart.err."//windowinnerdate//".nc"
            tmp_CO2nc_r4(:,:,:,1) = real(X_aErr3DCO2(:,:,:)/ppm )
            call Write_CO2_nc( filename_in, windowinnerdate, tmp_CO2nc_r4 )

            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine CO2_update


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Flx_extract(NowWin,X_a3DFlx,X_aErr3DFlx)
    !!=========================================================================
    !!
    !!  flux update this win and next win
    !!
    !!=========================================================================

        implicit none 
        !------------------------------------------------------------------

            !----------------------------------------------------------
            type(WindowsInfo),intent(in)    ::  NowWin
            real(8),intent(in)              ::  X_a3DFlx(nx,ny,nCoeRe)
            real(8),intent(in)              ::  X_aErr3DFlx(nx,ny,nCoeRe)

            !   loop
            integer                         ::  ii, jj, kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   read bpch
            !----------------------------------------------------------
            !   units: g C/m2/day                                     !
            !----------------------------------------------------------
            filename_in = trim(TTRunDIR)//"/Backrun/DigFlux.bpch"
            write(*,*) filename_in
            call Extract_Flux( Winlength, filename_in, tmp_Fluxwinlen_r4)

            !   output nc
            do jj=1,Winlength
                call day_forward(NowWin%WinStyear,NowWin%WinStmonth,NowWin%WinStday,&
                    &jj,WinInyear,WinInmonth,WinInday)
                write(windowinnerdate,'(I8)')  &
                    & WinInyear* 10000 + WinInmonth * 100 + WinInday
                !   Xb
                filename_in = trim(TTRunDIR)//"/Backrun"
                filename_in = trim(filename_in)//"/Flux."//windowinnerdate//".nc"
                tmp_Flux_r4(:,:,:,1) = tmp_Fluxwinlen_r4(:,:,:,jj)
                call Write_Flux_nc( filename_in, windowinnerdate, tmp_Flux_r4 )

                !   Xa
                filename_in = trim(TTRunDIR)//"/Backrun"
                filename_in = trim(filename_in)//"/Flux.xa."//windowinnerdate//".nc"
                do kk=1,nCoeRe
                tmp_Flux_r4(:,:,kk,1) = tmp_Fluxwinlen_r4(:,:,kk,jj) * X_a3DFlx(:,:,kk)
                enddo
                call Write_Flux_nc( filename_in, windowinnerdate, tmp_Flux_r4 )

                !   err
                filename_in = trim(TTRunDIR)//"/Backrun"
                filename_in = trim(filename_in)//"/Flux.err."//windowinnerdate//".nc"
                do kk=1,nCoeRe
                tmp_Flux_r4(:,:,kk,1) = tmp_Fluxwinlen_r4(:,:,kk,jj) * X_aErr3DFlx(:,:,kk)
                enddo
                call Write_Flux_nc( filename_in, windowinnerdate, tmp_Flux_r4 )
            enddo
            !----------------------------------------------------------
            return
        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine Flx_extract

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==================================================================================
end MODULE DAInOut
