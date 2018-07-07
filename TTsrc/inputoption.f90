MODULE InputOption
!!=========================================================================
!!
!!  Input time and options (not only options)
!!
!!  define type of input option to conform all the opt will be used.
!!  
!!      Type WindowsInfo
!!
!!      Subroutine LoadOptions
!!
!!  note : IOswitch part need to modify.
!!
!!  note : Samplescale, Winlength, Loc_uv, maxval_sat, beta
!!         need parameter calibration.
!!
!!    Record of revision:
!!           Date:          Programmer:          Description of change:
!!     ------------       --------------        ------------------------
!!     2017.10.29               Hanrui             NV2.0 original code
!!     2018.03.28               Hanrui             NV2.0 Final version
!!
!!=========================================================================
    implicit none
    !------------------------------------------------------------------

        !----------------------------------------------------------
        !   aassimilation parameters
        !sampling parameter
        integer,parameter   ::  alpha           =   8
        integer,parameter   ::  beta            =   8

        !landuse areas
        integer,parameter   ::  nlda_cot        =   23
        integer,parameter   ::  nlda_lai        =   120
        integer,parameter   ::  nlda_sea        =   31

        !inflation parameter
        real(8),parameter   ::  InflationFlx    =   1.05d0
        real(8),parameter   ::  InflationCO2    =   1.0d0

        !scale factor
        real(8),parameter   ::  SamplescaleCot  =   0.0d0   !countries
        real(8),parameter   ::  SamplescaleLai  =   0.1d0   !LAI
        real(8),parameter   ::  SamplescaleSea  =   0.1d0   !Ocean

        !parameters in localization
        real(8),parameter   ::  Loc_uv          =   3000.0d0
        !real(8),parameter   ::  Loc_Flux        =   900.0d0

        !lag window length (unit: week)
        integer,parameter   ::  lagWindow_CO2   =   0
        integer,parameter   ::  lagWindow_FLX   =   2

        !sampling fold days (unit: day)
        !   windows folding propulsions to give better results
        logical,parameter   ::  Fold_IOswitch   =   .True.
        integer,parameter   ::  Foldlen         =   0

        !   integer Winlength
        integer,parameter   ::  Winlength       =   3

        !   integer NRENS need to be settled before
        integer,parameter   ::  NRENS           =   50

        !----------------------------------------------------------
        !   IO options
        !flux region update run
        !   True: region average (GRG)
        !   False: gird update (GRD)
        logical,parameter   ::  Flxrun_RegionUp         = .False.
        !flux map opt
        !   True: use (TTT)
        !   False: don't use (FFF)
        logical,parameter   ::  Mapcot_IOswitch         = .True.
        logical,parameter   ::  Maplai_IOswitch         = .True.
        logical,parameter   ::  Mapsea_IOswitch         = .True.
        !flux update way for next window
        !   True: 1/3(F-1+Fnow+1.0)
        !   False: Flxrun_RegionUp + cot/lai/sea
        logical,parameter   ::  FluxNextWavg_IOswitch   = .True.
        !CO2 last day DA.(T/F)
        !   True: CO2 DA only last day. move last back to sp.
        !   False: no CO2 DA. move last back to sp.
        logical,parameter   ::  CO2ENDTimeDA_IOswitch   = .True.

        !----------------------------------------------------------
        !   OBSinfo
        !   Tansat_level
        !   OCO2 and GOSAT level
        integer,parameter   ::  Tansat_level    =   21
        integer,parameter   ::  OCOGST_level    =   20
        !   maxval of OCO2 or GST data in a single day.
        integer,parameter   ::  maxval_sat      =   20000

        !   Satellite OCO2 and GOSAT
        TYPE, PUBLIC :: Satellite_obs_bmc
            integer         ::  year!YYYYMMDDhh ( hour is ok )
            integer         ::  month
            integer         ::  day
            integer         ::  hour
            real            ::  lon
            real            ::  lat
            real            ::  xco2                    !XCO2 Value (Xco2)
            real            ::  xco2_err                !XCO2 uncertainty (err)
            integer         ::  xco2_wal                !Data Quality Indicator(warn level)
            real            ::  xco2_apriori            !A priori XCO2 Value (Xa)
            real            ::  xco2_avk(OCOGST_level)  !XCO2 Column Averaging Kernel (a)
            real            ::  co2_prof(OCOGST_level)  !CO2 Apriori Profile (Uap)
            real            ::  pres_wgh(OCOGST_level)  !Pressure Weighting Function (h)
            real            ::  pres_lev(OCOGST_level)  !Pressure Levels (P)
            integer         ::  reject                  !after thinning and quality Ctl
                                                        !   0 accept to use
                                                        !   1 reject to use
        END TYPE Satellite_obs_bmc

        !   Satellite TANSAT
        TYPE, PUBLIC :: TanSat_obs_bmc
            integer         ::  year!YYYYMMDDhh ( hour is ok )
            integer         ::  month
            integer         ::  day
            integer         ::  hour
            real            ::  lon
            real            ::  lat
            real            ::  xco2                    !XCO2 Value (Xco2)
            real            ::  xco2_err                !XCO2 uncertainty (err)
            real            ::  xco2_apriori            !A priori XCO2 Value (Xa)
            real            ::  xco2_avk(Tansat_level)  !XCO2 Column Averaging Kernel (a)
            real            ::  co2_prof(Tansat_level)  !CO2 Apriori Profile (Uap)
            real            ::  pres_wgh(Tansat_level)  !Pressure Weighting Function (h)
            real            ::  pres_lev(Tansat_level)  !Pressure Levels (P)
            integer         ::  reject                  !after thinning and quality Ctl
                                                        !   0 accept to use
                                                        !   1 reject to use
        END TYPE TanSat_obs_bmc

        !   OBSinfo after obs_operator
        TYPE, PUBLIC :: DA_obs_bmc
            real(8)         ::  lon
            real(8)         ::  lat
            real(8)         ::  alt
            real(8)         ::  obs                     !XCO2 Value (Xco2)
            real(8)         ::  obs_err                 !XCO2 uncertainty (err)
            real(8)         ::  obs_Xb
            real(8)         ::  obs_AA(NRENS)
        END TYPE DA_obs_bmc

        !----------------------------------------------------------
        !   windows information
        TYPE, PUBLIC :: WindowsInfo
            !windows time information
            integer         ::  WCYCLE
            integer         ::  WinStyear
            integer         ::  WinStmonth
            integer         ::  WinStday
            logical         ::  IFstart
            logical         ::  IFend
        END TYPE WindowsInfo

        !   inner time information
        integer             ::  WinInyear       !inner window year
        integer             ::  WinInmonth      !inner window month
        integer             ::  WinInday        !inner window day

        !   date information
        character(8)        ::  windowStdate
        character(8)        ::  windowEddate
        character(8)        ::  windowinnerdate

        !----------------------------------------------------------
        !   namelist
        !DIRCTORYs:
        character(80)       ::  TTDRIVERDIR
        character(80)       ::  TTRunDIR
        character(80)       ::  TTDATADIR
        !give NRENS to shell
        integer             ::  NRENStmp
        !Shell Cycle num:
        integer             ::  NOWCYCLE
        !NLS 5 Cycle num:
        integer             ::  NLS5num
        integer             ::  NLS5CYCLE
        !pre TIMEs and Parameters:
        integer             ::  STARTYEAR
        integer             ::  STARTMONTH
        integer             ::  ENDYEAR
        integer             ::  ENDMONTH
        !Samples and Forcing Ctl:
        character(8)        ::  StartTime
        character(8)        ::  EndTime
        integer             ::  Stflag
        integer             ::  Edflag

        namelist /TTnml/ &
            !DIRCTORYs:
            &TTDRIVERDIR,TTRunDIR,TTDATADIR,&
            !Shell Cycle num:
            &NRENStmp,NOWCYCLE,&
            !NLS 5 Cycle num:
            &NLS5num,NLS5CYCLE,&
            !pre TIMEs and Parameters:
            &STARTYEAR,STARTMONTH,ENDYEAR,ENDMONTH,&
            !Samples and Forcing time Ctl:
            &StartTime,EndTime,Stflag,Edflag
        !----------------------------------------------------------

    !------------------------------------------------------------------

    CONTAINS

    subroutine LoadOptions()
    !!=================================================================
    !!  assign directory , Time Ctrl and inputoptions
    !!=================================================================
        implicit none
        !----------------------------------------------------------

            !--------------------------------------------------
            integer     ::  unit
            unit = 55

            open(unit,file="TTnamelist")
            read(unit,nml=TTnml)
            !rewind(unit)    !start from the begining of the file
            close(unit)

            !--------------------------------------------------

        !----------------------------------------------------------

    !!=================================================================
    end subroutine LoadOptions

!!=========================================================================
end MODULE InputOption
