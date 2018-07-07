MODULE Common
!!=================================================================================
!!
!!    Input parameters
!!
!!    note: rx,ry,rz need parameter calibration.
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2017.10.29                   Hanrui                 NV2.0 original code
!!     2018.03.28                   Hanrui                 NV2.0 Final version
!!
!!=================================================================================
    IMPLICIT NONE

    PUBLIC
    !--------------------------------------------------------------------------
    
        !------------------------------------------------------------------
        ! Variable size definitions
        integer,parameter   ::  SHR_KIND_R8 = selected_real_kind(12) ! 8 byte
        integer,parameter   ::  SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte
        integer,parameter   ::  SHR_KIND_RN = kind(1.0)              ! native
        integer,parameter   ::  SHR_KIND_I8 = selected_int_kind (13) ! 8 byte
        integer,parameter   ::  SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte
        integer,parameter   ::  SHR_KIND_IN = kind(1)                ! native

        !------------------------------------------------------------------
   
    !--------------------------------------------------------------------------

        !------------------------------------------------------------------
        !model girds
        integer,parameter       ::  nx              =   144
        integer,parameter       ::  ny              =   91
        integer,parameter       ::  nz              =   47
        integer,parameter       ::  nCoeRe          =   10
        !loc mod
        integer,parameter       ::  rx              =   30
        integer,parameter       ::  ry              =   30
        integer,parameter       ::  rz              =   1 ! dont do loc on zdir
        ! one day
        integer,parameter       ::  nt              =   1
        ! nbox is 2*2 box to itp nbox must eq to2
        integer,parameter       ::  nbox            =   2
        ! ntimes is times of DigInst var (co2,pres) in one day.
        integer,parameter       ::  ntimes          =   8

        !Natural parameters
        real(8),parameter       ::  ppm             =   1.0E6
        real(8),parameter       ::  na              =   6.0225d0*1.0E23
        real(8),parameter       ::  PI              =   3.1415926535898d0
        real(8),parameter       ::  Rearth          =   6371.004d0
        real(8),parameter       ::  km              =   1000.0d0
        real(8),parameter       ::  sec_in_year1    =   86400.0d0*365.0d0

        !file io parameter (subroutine give_empty_unit/ntp)
        integer,parameter       ::  unit1           =   11
        integer,parameter       ::  unit2           =   22
        integer,parameter       ::  unit3           =   33
        integer,parameter       ::  unit4           =   44

        !------------------------------------------------------------------
        !do cycle parameter 
        !   bus error
        !   when sub1() call sub2() in do ii cycle, and sub2() use ii.
        ! PRIVATE
        integer,PRIVATE         ::  i,      j,      k,      l 
        integer,PRIVATE         ::  ii,     jj,     kk,     ll 

        !character parameters
        integer,parameter       ::  Monthdays(12)   =   &
                    &(/31,28,31,30,31,30,31,31,30,31,30,31/)

        character*2,parameter   ::  NumCha(100)     =   &
                    &(/'01','02','03','04','05','06','07','08','09','10',&
                      &'11','12','13','14','15','16','17','18','19','20',&
                      &'21','22','23','24','25','26','27','28','29','30',&
                      &'31','32','33','34','35','36','37','38','39','40',&
                      &'41','42','43','44','45','46','47','48','49','50',&
                      &'51','52','53','54','55','56','57','58','59','60',&
                      &'61','62','63','64','65','66','67','68','69','70',&
                      &'71','72','73','74','75','76','77','78','79','80',&
                      &'81','82','83','84','85','86','87','88','89','90',&
                      &'91','92','93','94','95','96','97','98','99','00'/)

        !   character date YYYYMMDD
        character(8)            ::  DateCha

        !gird information
        !   longitude
        real(8),parameter       ::  MLonC(nx)       =   &
                &(/-180.0d0,-177.5d0,-175.0d0,-172.5d0,-170.0d0,-167.5d0, &
                  &-165.0d0,-162.5d0,-160.0d0,-157.5d0,-155.0d0,-152.5d0, &
                  &-150.0d0,-147.5d0,-145.0d0,-142.5d0,-140.0d0,-137.5d0, &
                  &-135.0d0,-132.5d0,-130.0d0,-127.5d0,-125.0d0,-122.5d0, &
                  &-120.0d0,-117.5d0,-115.0d0,-112.5d0,-110.0d0,-107.5d0, &
                  &-105.0d0,-102.5d0,-100.0d0, -97.5d0, -95.0d0, -92.5d0, &
                  & -90.0d0, -87.5d0, -85.0d0, -82.5d0, -80.0d0, -77.5d0, &
                  & -75.0d0, -72.5d0, -70.0d0, -67.5d0, -65.0d0, -62.5d0, &
                  & -60.0d0, -57.5d0, -55.0d0, -52.5d0, -50.0d0, -47.5d0, &
                  & -45.0d0, -42.5d0, -40.0d0, -37.5d0, -35.0d0, -32.5d0, &
                  & -30.0d0, -27.5d0, -25.0d0, -22.5d0, -20.0d0, -17.5d0, &
                  & -15.0d0, -12.5d0, -10.0d0,  -7.5d0,  -5.0d0,  -2.5d0, &
                  &   0.0d0,   2.5d0,   5.0d0,   7.5d0,  10.0d0,  12.5d0, &
                  &  15.0d0,  17.5d0,  20.0d0,  22.5d0,  25.0d0,  27.5d0, &
                  &  30.0d0,  32.5d0,  35.0d0,  37.5d0,  40.0d0,  42.5d0, &
                  &  45.0d0,  47.5d0,  50.0d0,  52.5d0,  55.0d0,  57.5d0, &
                  &  60.0d0,  62.5d0,  65.0d0,  67.5d0,  70.0d0,  72.5d0, &
                  &  75.0d0,  77.5d0,  80.0d0,  82.5d0,  85.0d0,  87.5d0, &
                  &  90.0d0,  92.5d0,  95.0d0,  97.5d0, 100.0d0, 102.5d0, &
                  & 105.0d0, 107.5d0, 110.0d0, 112.5d0, 115.0d0, 117.5d0, &
                  & 120.0d0, 122.5d0, 125.0d0, 127.5d0, 130.0d0, 132.5d0, &
                  & 135.0d0, 137.5d0, 140.0d0, 142.5d0, 145.0d0, 147.5d0, &
                  & 150.0d0, 152.5d0, 155.0d0, 157.5d0, 160.0d0, 162.5d0, &
                  & 165.0d0, 167.5d0, 170.0d0, 172.5d0, 175.0d0, 177.5d0/) 

        !   latitude
        real(8),parameter       ::  MLatC(ny)       =   &
                &(/-89.5d0,-88.0d0,-86.0d0,-84.0d0,-82.0d0,-80.0d0,-78.0d0,&
                  &-76.0d0,-74.0d0,-72.0d0,-70.0d0,-68.0d0,-66.0d0,-64.0d0,&
                  &-62.0d0,-60.0d0,-58.0d0,-56.0d0,-54.0d0,-52.0d0,-50.0d0,&
                  &-48.0d0,-46.0d0,-44.0d0,-42.0d0,-40.0d0,-38.0d0,-36.0d0,&
                  &-34.0d0,-32.0d0,-30.0d0,-28.0d0,-26.0d0,-24.0d0,-22.0d0,&
                  &-20.0d0,-18.0d0,-16.0d0,-14.0d0,-12.0d0,-10.0d0, -8.0d0,&
                  & -6.0d0, -4.0d0, -2.0d0,  0.0d0,  2.0d0,  4.0d0,  6.0d0,&
                  &  8.0d0, 10.0d0, 12.0d0, 14.0d0, 16.0d0, 18.0d0, 20.0d0,&
                  & 22.0d0, 24.0d0, 26.0d0, 28.0d0, 30.0d0, 32.0d0, 34.0d0,&
                  & 36.0d0, 38.0d0, 40.0d0, 42.0d0, 44.0d0, 46.0d0, 48.0d0,&
                  & 50.0d0, 52.0d0, 54.0d0, 56.0d0, 58.0d0, 60.0d0, 62.0d0,&
                  & 64.0d0, 66.0d0, 68.0d0, 70.0d0, 72.0d0, 74.0d0, 76.0d0,&
                  & 78.0d0, 80.0d0, 82.0d0, 84.0d0, 86.0d0, 88.0d0, 89.5d0/)

        !   level num (use in creating nc)
        real(8),parameter       ::  MLevC(nz)       =   &
                &   (/ (kk*1.0d0, kk = 1, nz) /)

        !   Flux level num (use in creating nc)
        real(8),parameter       ::  MCoeReC(nCoeRe) =   &
                &   (/ (ll*1.0d0, ll = 1, nCoeRe) /)

        !   Flux level name (use as name in creating nc)
        character(len=6)        ::  FluxSpi(nCoeRe) =   &
                &   (/ "CO2ff " ,"CO2oc " ,"CO2bal" ,"CO2bb " ,"CO2bf " ,&
                &      "CO2nte" ,"CO2shp" ,"CO2pln" ,"CO2che" ,"CO2sur" /)

        !   altitude km (obsoleted, as a referance info)
        real(8),parameter       ::  LevAlt(nz)      =   &
                &(/ 0.058d0, 0.189d0, 0.320d0, 0.454d0, 0.589d0, 0.726d0,&
                  & 0.864d0, 1.004d0, 1.146d0, 1.290d0, 1.436d0, 1.584d0,&
                  & 1.759d0, 1.988d0, 2.249d0, 2.517d0, 2.792d0, 3.074d0,&
                  & 3.439d0, 3.896d0, 4.375d0, 4.879d0, 5.413d0, 5.980d0,&
                  & 6.585d0, 7.237d0, 7.943d0, 8.846d0, 9.936d0,11.021d0,&
                  &12.086d0,13.134d0,14.170d0,15.198d0,16.222d0,17.243d0,&
                  &18.727d0,20.836d0,23.020d0,25.307d0,28.654d0,34.024d0,&
                  &40.166d0,47.135d0,54.834d0,63.058d0,72.180d0/)   !km

        !   altitude hPa (obsoleted)
        real(8),parameter       ::  LevPre(nz)      =   &
                &(/1005.0d0, 990.0d0, 975.0d0, 959.0d0, 944.0d0, 929.0d0,&
                  & 913.0d0, 898.0d0, 883.0d0, 868.0d0, 852.0d0, 837.0d0,&
                  & 819.0d0, 796.0d0, 771.0d0, 745.0d0, 720.0d0, 694.0d0,&
                  & 663.0d0, 624.0d0, 586.0d0, 548.0d0, 510.0d0, 472.0d0,&
                  & 434.0d0, 396.0d0, 358.0d0, 313.0d0, 267.0d0, 226.0d0,&
                  & 192.0d0, 163.0d0, 139.0d0, 118.0d0, 100.0d0,  85.0d0,&
                  &  67.0d0,  48.0d0,  34.0d0,  24.0d0,  14.0d0,   6.0d0,&
                  &   2.0d0,   1.0d0,   0.4d0,   0.1d0,  0.04d0/)   !hPa


    !--------------------------------------------------------------------------

!==================================================================================
end MODULE Common
