MODULE extract_bpch
!!=================================================================================
!!
!!  Description:
!!      Extract instantaneous CO2 and daily Flux from Model results.
!!      files:
!!          DigInst.YYYYMMDD.bpch
!!          DigFlux.bpch
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2018.01.02                   Hanrui                  NV2.0 original code
!!     2018.03.28                   Hanrui                  NV2.0 Final version
!!
!!=================================================================================

    !USE and IMPLICIT
    use Common,             only    :   nx,ny,nz,nCoeRe,ntimes
    use Common,             only    :   na,sec_in_year1,NumCha
    use inputoption,        only    :   Winlength
    implicit none

    !--------------------------------------------------------------------------
    !INTERFACE
    PUBLIC                  ::  Extract_InstCO2

    PUBLIC                  ::  Extract_Flux

    !read from r4 to r8
    real(4),PUBLIC          ::  tmp_CO2_r4(nx,ny,nz,ntimes)
    real(4),PUBLIC          ::  tmp_Prue_r4(nx,ny,nz,ntimes)
    real(4),PUBLIC          ::  tmp_Fluxwinlen_r4(nx,ny,nCoeRe,Winlength)

    !--------------------------------------------------------------------------
    !PRIVATE VARIABLES
    PRIVATE
        !------------------------------------------------------------------
        !inner private Variables for bpch file
        real(4)             ::  temp(nx,ny,nz+1,1)    
        real(4)             ::  gridarea(nx,ny,1,1) 
      
        character(40)       ::  fti
        character(80)       ::  title

        integer             ::  NNL,ios
        integer             ::  tCO2
        integer             ::  tPrue
        integer             ::  tFlx(10)
        integer             ::  tOth
        integer             ::  unit1

        integer             ::  NTRACER,   NSKIP
        integer             ::  NTID
        integer             ::  HALFPOLAR, CENTER180
        integer             ::  NI,        NJ,        NL
        integer             ::  IFIRST,    JFIRST,    LFIRST
        real(4)             ::  LONRES,    LATRES
        real(8)             ::  ZTAU0,     ZTAU1
        character(20)       ::  MODELNAME
        character(40)       ::  CATEGORY
        character(40)       ::  UNIT  
        character(40)       ::  RESERVED

        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    Contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Extract_InstCO2( InFile, CO2C, Prue)
    !!=========================================================================
    !!
    !!  Description: 
    !!      Extract CO2 Concentration from DigInst.YYYYMMDD.bpch
    !!      Daily instantaneous CO2, 8 times.
    !!      03 06 09 12 15 18 21 24
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in
            character(*),   intent(in)  ::  InFile
            !   out 
            real(4),        intent(out) ::  CO2C( nx, ny, nz, 8)
            real(4),        intent(out) ::  Prue( nx, ny, nz, 8)

            !   inner
            !   edge presure level = nz + 1
            real(4)                     ::  PrueEdge( nx, ny, nz+1, 8)

            !   loop
            integer             ::  i,  j,  k,  l
            integer             ::  ii, jj, kk, ll, tt

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   Set initial values
            unit1 = 11
            temp  = 0.0
            CO2C  = 0.0
            PrueEdge = 0.0
            tCO2  = 0
            tPrue = 0
            write(*,*) trim(infile)

            !----------------------------------------------------------
            !   open input files and output files
            open(unit1,file=trim(infile),convert='big_endian',&
                  &status='old',iostat=ios,form='unformatted')
            if( ios /= 0 )                 print*,'OPEN FILE ERROR'

            !----------------------------------------------------------
            !   read_header
            read(unit1,iostat=ios)fti
            if( ios /= 0 )                 print*,'READ RTI ERROR'
            if(trim(fti)/='CTM bin 02')    print*,'INVALID FILE FORMAT'
            if( ios == 0 )                 print*,fti
            read(unit1,iostat=ios)title
            if( ios /= 0 )                 print*,'OPEN TITLE ERROR'
            if( ios == 0 )                 print*,title   

        do

            !----------------------------------------------------------
            !   Read Var                                              !
            !----------------------------------------------------------
            read( unit1,iostat=ios) MODELNAME, LONRES, LATRES,&
                & HALFPOLAR, CENTER180
            if   ( ios < 0 )         exit
            if   ( ios > 0 )         print*,'READ MODELNAME ERROR'
            read( unit1,iostat=ios) CATEGORY, NTRACER, UNIT, ZTAU0, ZTAU1, &
            &RESERVED, NI, NJ, NL,IFIRST, JFIRST, LFIRST, NSKIP
            if   ( ios /= 0)         print*,'READ CTEGORY ERROR'
            read(unit1,iostat=ios) (((temp(i,j,k,1),i=1,NI),j=1,NJ),k=1,NL)

            !----------------------------------------------------------
            !   Extracting CO2 concentration                          !
            !----------------------------------------------------------
            if(trim(CATEGORY)==trim('IJ-AVG-$').and.NTRACER==1)then
            !   Time loop
            tCO2 = tCO2 + 1
            if( LFIRST/=1 )        print*,'LFIRST/=1-11'
            NNL=NL+LFIRST-1
                  CO2C( 1:NI, 1:NJ, 1:NZ, tCO2:tCO2) = &
                & temp( 1:NI, 1:NJ, 1:NZ, 1:1)
            endif

            !----------------------------------------------------------
            !   Extracting gird edge Presure                          !
            !----------------------------------------------------------
            if(trim(CATEGORY)==trim('PEDGE-$').and.NTRACER==1)then
            !   Time loop
            tPrue = tPrue + 1
            if( LFIRST/=1 )        print*,'LFIRST/=1-11'
            NNL=NL+LFIRST-1
                  PrueEdge( 1:NI, 1:NJ, 1:NZ+1, tPrue:tPrue) = &
                &     temp( 1:NI, 1:NJ, 1:NZ+1, 1:1)
            endif

        enddo

            !----------------------------------------------------------
            !   clean
            close(unit1)

            !----------------------------------------------------------
            !   Calculate Prue
            do kk = 1, nz
                Prue(:,:,kk,:) = PrueEdge(:,:,kk,:) + PrueEdge(:,:,kk+1,:)
                Prue(:,:,kk,:) = Prue(:,:,kk,:) / 2.0
            end do

            !----------------------------------------------------------
            return

        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine Extract_InstCO2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Extract_Flux( nt, InFile, Flux)
    !!=======================================================================80
    !
    !   extract flux from GEOS-Chem OUTPUT  :
    !
    !   VAR      NAME                           NTRACER         NLambda     
    !   CO2ff    CO2 fossil fuel emiss                1               1 (cot)
    !   CO2oc    CO2 ocean emissions                  2               3 (sea)
    !   CO2bal   CO2 balanced biosphere               3               2 (lai)
    !   CO2bb    CO2 biomass burning emiss            4           4 NaN .
    !   CO2bf    CO2 biofuel emission                 5           4 NaN .
    !   CO2nte   CO2 net terr exchange                6               2 (lai)
    !   CO2shp   CO2 ship emissions                   7           4 NaN .
    !   CO2pln   CO2 aircraft emissions               8(3D)       4 NaN .
    !   CO2che   CO2 chemical oxidation               9(3D)       4 NaN .
    !   CO2sur   CO2 chem surf correction            10           4 NaN (oth)
    !
    !!=======================================================================80

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in
            integer,        intent(in)  ::  nt
            character(*),   intent(in)  ::  InFile
            !   out 
            real(4),        intent(out) ::  Flux( nx, ny, nCoeRe, nt)

            !----------------------------------------------------------
            !   Calculate na=6.0225E23 must use r8
            real(8)                     ::  Flux_r8( nx, ny, nCoeRe)

            !   loop
            integer             ::  i,  j,  k,  l
            integer             ::  ii, jj, kk, ll, tt

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   Set initial values
            unit1 = 11
            Flux  = 0.0
            temp  = 0.0
            gridarea = 0.0
            tFlx  = 0
            tOth  = 0
            write(*,*) trim(InFile)

            !----------------------------------------------------------
            !   open input files and output files
            open(unit1,file=trim(InFile),convert='big_endian',&
                  &status='old',iostat=ios,form='unformatted')
            if( ios /= 0 )                 print*,'OPEN FILE ERROR'

            !----------------------------------------------------------
            !   read_header
            read(unit1,iostat=ios)fti
            if( ios /= 0 )                 print*,'READ RTI ERROR'
            if(trim(fti)/='CTM bin 02')    print*,'INVALID FILE FORMAT'
            if( ios == 0 )                 print*,fti

            read(unit1,iostat=ios)title
            if( ios /= 0 )                 print*,'OPEN TITLE ERROR'
            if( ios == 0 )                 print*,title    

        do

            !----------------------------------------------------------
            !   Read Var                                              !
            !----------------------------------------------------------
            read( unit1,iostat=ios) MODELNAME, LONRES, LATRES,&
                & HALFPOLAR, CENTER180
            if   ( ios < 0 )         exit
            if   ( ios > 0 )         print*,'READ MODELNAME ERROR'
            read( unit1,iostat=ios) CATEGORY, NTRACER, UNIT, ZTAU0, ZTAU1, &
            &RESERVED, NI, NJ, NL,IFIRST, JFIRST, LFIRST, NSKIP
            if   ( ios /= 0)         print*,'READ CTEGORY ERROR'
            read(unit1,iostat=ios) (((temp(i,j,k,1),i=1,NI),j=1,NJ),k=1,NL)

            !----------------------------------------------------------
            !   Extracting Flux (atmos C/cm2/s)                       !
            !----------------------------------------------------------
            do NTID = 1, nCoeRe

            if(trim(CATEGORY) == trim('CO2-SRCE') .and. NTRACER == NTID) then
            !   Time loop
            tFlx(NTID) = tFlx(NTID) + 1

            if( LFIRST/=1 )  print*,'LFIRST/=1-'//NumCha(NTID)
            NNL=NL+LFIRST-1

            if(NTRACER == 8 .or. NTRACER ==9) then
                do kk = 1,nz
                  Flux( 1:NI, 1:NJ, NTID:NTID, tFlx(NTID):tFlx(NTID)) = &
                & Flux( 1:NI, 1:NJ, NTID:NTID, tFlx(NTID):tFlx(NTID)) + &
                & temp( 1:NI, 1:NJ, kk:kk,     1:1)
                end do
            else
                  Flux( 1:NI, 1:NJ, NTID:NTID, tFlx(NTID):tFlx(NTID)) = &
                & temp( 1:NI, 1:NJ, 1:1,       1:1)
            endif

            endif

            enddo

            !----------------------------------------------------------
            !   Extrating area m2                                     !
            !----------------------------------------------------------
            if(tOth == 0) then
            if(trim(CATEGORY)==trim('DXYP').and.NTRACER==1)then
            if( LFIRST/=1 )        print*,'LFIRST/=1-11'
            NNL=NL+LFIRST-1
            tOth = tOth + 1
              gridarea( 1:NI, 1:NJ, 1:1, tOth:tOth) = &
                & temp( 1:NI, 1:NJ, 1:1, 1:1)
            endif
            endif

        end do

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   Calculate (atmos C/cm2/s-->g C/m2/day)                !
            !----------------------------------------------------------
            do tt = 1, nt

            Flux_r8 = 0.0d0

            do kk = 1, nCoeRe
            do jj = 1, ny
            do ii = 1, nx

            Flux_r8( ii, jj, kk) = Flux( ii, jj, kk, tt) * &
                        & 1.0E4 *        &              !/cm2    -> /m2
                        & 60 * 60 * 24 * &              !s-1     -> day-1
                        & 12*1.0d0 /na                  !atmos C -> g
                        !& gridarea( ii, jj, 1, 1) * &   !/m2     -> /gird
                        !& 1.0E-3 *       &              !g       -> kg

            end do
            end do
            end do

            Flux(:,:,:,tt) = real( Flux_r8(:,:,:), 4)

            end do

            !----------------------------------------------------------
            return

        !------------------------------------------------------------------

    !!=========================================================================
    end subroutine Extract_Flux

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==================================================================================
end module extract_bpch
