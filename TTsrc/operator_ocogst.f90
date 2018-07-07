MODULE operator_ocogst
!!=================================================================================
!!
!!  Description:
!!      OBS operator of oco2.
!!
!!  Warning:
!!      var(lon,lat,lev) in nc_file, into var(lev,lat,lon) in fortran.
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2018.01.02                   Hanrui                  NV2.0 original code
!!     2018.03.28                   Hanrui                  NV2.0 Final version
!!
!!=================================================================================

    use common
    use netcdf

    use InputOption,    only :  NRENS
    use InputOption,    only :  OCOGST_level
    use InputOption,    only :  maxval_sat
    use InputOption,    only :  Satellite_obs_bmc
    use InputOption,    only :  DA_obs_bmc

    use interpolation,  only :  plane_4lagr,line_4lagr
    use interpolation,  only :  plane_3lagr,line_3lagr
    use interpolation,  only :  plane_2lagr,line_2lagr

    use extract_bpch,   only :  Extract_InstCO2
    use extract_bpch,   only :  tmp_CO2_r4,tmp_Prue_r4

    use Statistic,      only :  XSigma

    implicit none
    include "netcdf.inc"
    !--------------------------------------------------------------------------
    !   find OCO2 or gosat file name and length
    PUBLIC                  ::  OCOGST_Filename_length
    !   read in OCO2 or gosat obs
    PUBLIC                  ::  OCOGST_read
    !   OCO2 or gosat obs quality ctl
    PUBLIC                  ::  OCOGST_QualityCtl
    !   OCO2 or gosat obs operator
    PUBLIC                  ::  OCOGST_Operator
    !   OCO2 or gosat obs operator from CO2 and pres to Simulated XCO2(:)
    PUBLIC                  ::  OCOGST_Optcore
    !   locate start point of long and lat
    PUBLIC                  ::  locate_startpoint
    PUBLIC                  ::  locate_Boxvalue
    PUBLIC                  ::  locate_Boxaxis
    PUBLIC                  ::  locate_weight
    !   handle error
    PUBLIC                  ::  handle_err

    !--------------------------------------------------------------------------
    !PRIVATE VARIABLES
    PRIVATE

        !------------------------------------------------------------------
        !inner var:
        !   Attributes char
        character(len=128)      ::  AName
        !   file id and Status
        integer                 ::  Status
        integer                 ::  ncid
        !   variable id 
        integer                 ::  varid
        integer                 ::  dimid

        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    Contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine OCOGST_Filename_length( obsdate, filename_in, length, AChar)
    !==========================================================================
    !!
    !!  find OCO2 file name and length
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(8),intent(in)         ::  obsdate
        character(*),intent(inout)      ::  filename_in
        integer,intent(out)             ::  length
        character(3),intent(in)         ::  AChar


        character(120)                  ::  filename
        character(120)                  ::  namechar
        integer                         ::  Status
        !   loop
        integer                         ::  ii,     jj,     kk

        !------------------------------------------------------------------
        if(AChar.eq."oco") then
            filename_in = trim(filename_in)//"/Obs/OCO2_L2_Lite_FP.8r/"
        elseif(AChar.eq."gst") then
            filename_in = trim(filename_in)//"/Obs/ACOS_L2_Lite_FP.7.3/"
        else
            stop "only oco and gst is allowed here in OCOGST_Filename_length"
        endif
        filename = trim(filename_in)//obsdate(1:4)//"/"//"name.txt"

        !   file exist          length = length
        !   file don't exist    length = -1 
        length = 1

        open(unit1,file=trim(filename),iostat=Status)
        if ( Status /= 0 ) then
            write(*,*) "Error opening file name :",trim(filename)
            length = -1
        endif
        
        do while( length.eq.1 )

            read(unit1,*,iostat=Status) namechar
            !write(*,*) namechar
            !write(*,*) namechar(12:17),obsdate(3:8)

            if(Status.lt.0) then
                length=-1
                exit
            endif

            if(namechar(12:17).eq.obsdate(3:8)) then
                filename_in = trim(filename_in)//obsdate(1:4)//"/"//trim(namechar)
                exit
            endif

        enddo
        close(unit1)

        if(length.eq.1) then

            !OPen NETCDF file
            Status = nf_open( trim(filename_in), NF_NOWRITE, ncid)
            if (Status /= NF_NOERR) then
            write(*,*) "OPEN ERROR CHECK_FILE: ", trim(filename_in)
            stop
            end if

            !inquire varid dimension
            AName  = "sounding_id"
            Status = nf_inq_varid( ncid, trim(AName), varid)
            call handle_err( Status, filename_in, AName)

            !read var
            Status = nf_inq_dimlen( ncid, varid, length)
            call handle_err( Status, filename_in, AName)

            !close
            Status = nf_close( ncid)

            !write(*,*) trim(filename_in)
            !write(*,*) length

        endif

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine OCOGST_Filename_length

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine OCOGST_read( filename_in, length, daily_OCO2)
    !==========================================================================
    !!
    !!  read in OCO2 obs
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(*),intent(in)                 ::  filename_in
        integer,intent(in)                      ::  length
        type(Satellite_obs_bmc),intent(out)     ::  daily_OCO2(length)

        !------------------------------------------------------------------
        !   time
        integer                         ::  epoch
        integer(2),allocatable          ::  tmp_int2(:,:)
        integer(1)                      ::  tmp_int1(length)

        !   lon,lat,
        real                            ::  tmp_flot(length)
        real                            ::  tmp_flt2(OCOGST_level,length)
        !   loop
        integer                         ::  ii,     jj,     kk

        !------------------------------------------------------------------

        !OPen NETCDF file
        write(*,*) trim(filename_in)
        Status = nf_open( trim(filename_in), NF_NOWRITE, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) "OPEN ERROR CHECK_FILE: ", trim(filename_in)
        stop
        end if

        !read epoch dimension length
        AName  = "epoch_dimension"
        Status = nf_inq_dimid( ncid, trim(AName), dimid)
        call handle_err( Status, filename_in, AName)
        Status = nf_inq_dimlen( ncid, dimid, epoch)
        call handle_err( Status, filename_in, AName)
        !write(*,*) " epoch",epoch
        !write(*,*) "length",length

        !read time
        allocate(tmp_int2(epoch,length))
        AName  = "date"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_int2( ncid, varid, tmp_int2)
        call handle_err( Status, filename_in, AName)
        ! write(*,*) tmp_int2(:,10000)
        ! write(*,*) tmp_int2(:,20000)
        ! write(*,*) tmp_int2(:,length)
        do ii=1,length
            daily_OCO2(ii)%year     = tmp_int2(1,ii)
            daily_OCO2(ii)%month    = tmp_int2(2,ii)
            daily_OCO2(ii)%day      = tmp_int2(3,ii)
            daily_OCO2(ii)%hour     = tmp_int2(4,ii)
            ! hour from 00 to 23, hour = hour + 1 is not better
            ! write(*,*) daily_OCO2(ii)%year,daily_OCO2(ii)%month,&
            ! &daily_OCO2(ii)%day,daily_OCO2(ii)%hour
        enddo
        
        !read longitude
        AName  = "longitude"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flot(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%lon      = tmp_flot(ii)
            ! write(*,*) daily_OCO2(ii)%lon
        enddo
        
        !read latitude
        AName  = "latitude"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flot(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%lat      = tmp_flot(ii)
            ! write(*,*) daily_OCO2(ii)%lat
        enddo

        !read xco2
        AName  = "xco2"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flot(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%xco2     = tmp_flot(ii)
            ! write(*,*) daily_OCO2(ii)%xco2
        enddo
        
        !read xco2_err
        AName  = "xco2_uncertainty"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flot(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%xco2_err = tmp_flot(ii)
            ! write(*,*) daily_OCO2(ii)%xco2_err
        enddo
        
        !read warn_level
        AName  = "warn_level"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_int1(ncid, varid, tmp_int1(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%xco2_wal = tmp_int1(ii)
            ! write(*,*) daily_OCO2(ii)%xco2_wal
        enddo
        
        !read xco2_apriori
        AName  = "xco2_apriori"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flot(:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%xco2_apriori = tmp_flot(ii)
            ! write(*,*) daily_OCO2(ii)%xco2_apriori
        enddo

        !read xco2_averaging_kernel
        AName  = "xco2_averaging_kernel"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flt2(:,:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%xco2_avk(:) = tmp_flt2(:,ii)
            ! write(*,*) daily_OCO2(ii)%xco2_avk(:)
        enddo

        !read co2_profile_apriori
        AName  = "co2_profile_apriori"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flt2(:,:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%co2_prof(:) = tmp_flt2(:,ii)
            ! write(*,*) daily_OCO2(ii)%co2_prof(:)
        enddo

        !read pressure_weight
        AName  = "pressure_weight"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flt2(:,:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%pres_wgh(:) = tmp_flt2(:,ii)
            ! write(*,*) daily_OCO2(ii)%pres_wgh(:)
        enddo

        !read pressure_levels
        AName  = "pressure_levels"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_real(ncid, varid, tmp_flt2(:,:))
        call handle_err( Status, filename_in, AName)
        do ii=1,length
            daily_OCO2(ii)%pres_lev(:) = tmp_flt2(:,ii)
            ! write(*,*) daily_OCO2(ii)%pres_lev(:)
        enddo

        !close
        Status = nf_close( ncid)


        !check
        ! do jj=1,10
        !     write(*,*) "year:",daily_OCO2(jj)%year
        !     write(*,*) "month:",daily_OCO2(jj)%month
        !     write(*,*) "day:",daily_OCO2(jj)%day
        !     write(*,*) "hour:",daily_OCO2(jj)%hour
        !     write(*,*) "lon:",daily_OCO2(jj)%lon
        !     write(*,*) "lat:",daily_OCO2(jj)%lat
        !     write(*,*) "xco2:",daily_OCO2(jj)%xco2
        !     write(*,*) "xco2_err:",daily_OCO2(jj)%xco2_err
        !     write(*,*) "xco2_wal:",daily_OCO2(jj)%xco2_wal
        !     write(*,*) "xco2_apriori:",daily_OCO2(jj)%xco2_apriori
        !     write(*,*) "xco2_avk:",daily_OCO2(jj)%xco2_avk
        !     write(*,*) "co2_prof:",daily_OCO2(jj)%co2_prof
        !     write(*,*) "pres_wgh:",daily_OCO2(jj)%pres_wgh
        !     write(*,*) "pres_lev:",daily_OCO2(jj)%pres_lev
        !     write(*,*) "reject:",daily_OCO2(jj)%reject
        !     pause
        ! enddo
        
        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine OCOGST_read

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine OCOGST_QualityCtl( obsdate, filename_in, length, daily_OCO2, AChar)
    !==========================================================================
    !!
    !!  OCO2 obs quality ctl
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(8),intent(in)                 ::  obsdate
        character(*),intent(in)                 ::  filename_in
        !   filename_in pwd of TTRun
        integer,intent(in)                      ::  length
        type(Satellite_obs_bmc),intent(inout)   ::  daily_OCO2(length)
        character(3),intent(in)                 ::  AChar
        !   loop
        integer             ::  ii,     jj,     kk

        integer             ::  totalnum
        integer             ::  leap

        !   inner var
        !   back results
        real(8)             ::  Bgresult(length)
        integer             ::  Backlen
        real(8),allocatable ::  Backres(:)
        real(8)             ::  STDBack
        character(120)      ::  filename
        real(8)             ::  CO2_inst( nx, ny, nz, ntimes)
        real(8)             ::  PRES_inst( nx, ny, nz, ntimes)

        !------------------------------------------------------------------

        !   reject = 0 : accept
        !   reject = 1 : reject
        daily_OCO2(:)%reject = 0
        totalnum = length - sum(daily_OCO2(:)%reject)
        write(*,*) "Total num in one day:             ",totalnum

        !   data quality ctl
        !
        !   1 : warn_level 0(50%) 1(60%) 2(70%) 3(80%) 4(90%) 5(100%)
        !                  ^                           ^
        !                 OCO2                        GOSAT
        !
        !   2 : location ( range -88.0 ~ 88.0 )
        !
        !   3 : O-B test
        !
        !   4 : missing value (xco2_wal filtered)
        !
        if(AChar.eq."oco") then
            do ii=1,length
            if( daily_OCO2(ii)%xco2_wal.ne.0 .or. daily_OCO2(ii)%xco2_err.gt.2.0 ) then
                daily_OCO2(ii)%reject = 1
            endif
            if( daily_OCO2(ii)%lat.gt.88.0 .or. daily_OCO2(ii)%lat.lt.-88.0) then
                daily_OCO2(ii)%reject = 1
            endif
            enddo
        elseif(AChar.eq."gst") then
            do ii=1,length
            if( daily_OCO2(ii)%xco2_wal.gt.4 .or. daily_OCO2(ii)%xco2_err.gt.5.0 ) then
                daily_OCO2(ii)%reject = 1
            endif
            if( daily_OCO2(ii)%lat.gt.88.0 .or. daily_OCO2(ii)%lat.lt.-88.0 ) then
                daily_OCO2(ii)%reject = 1
            endif
            enddo
        else
            stop "only oco and gst is allowed here in OCOGST_QualityCtl"
        endif

        !   O-B test
        !   read in xb co2, and pres
        filename = trim(filename_in)//"/DigInst."//obsdate//".bpch"
        call Extract_InstCO2( filename, tmp_CO2_r4, tmp_Prue_r4)
        CO2_inst(:,:,:,:)   = tmp_CO2_r4(:,:,:,:) * 1.0d0 * ppm
        PRES_inst(:,:,:,:)  = tmp_Prue_r4(:,:,:,:) * 1.0d0

        call OCOGST_Optcore( CO2_inst, PRES_inst, length, daily_OCO2, Bgresult)

        Backlen = length - sum(daily_OCO2(:)%reject)
        allocate(Backres(Backlen))
        jj=0
        do ii=1,length
            if( daily_OCO2(ii)%reject .ne. 1) then
            jj = jj+1
            Backres(jj) = daily_OCO2(ii)%xco2 - Bgresult(ii)
            endif
        enddo
        STDBack = 3 * XSigma(Backlen,Backres)
        write(*,*) "O-B 3 STD err =",STDBack,"ppm"
        do ii=1,length
        !   empericial for real data assimilation only
        if( Bgresult(ii) - daily_OCO2(ii)%xco2 .gt. STDBack .and.&
          & Bgresult(ii) - daily_OCO2(ii)%xco2 .gt. 6.0d0 ) then
            daily_OCO2(ii)%reject = 1 
        endif
        enddo
        deallocate(Backres)

        !   data thinning
        totalnum = length - sum(daily_OCO2(:)%reject)
        write(*,*) "num in one day after quality ctl: ",totalnum

        if(totalnum.gt.maxval_sat) then

            leap = ceiling((totalnum*1.0)/(maxval_sat*1.0))
            !   ceiling : int upword
            write(*,*) "thinning leap use:                ",leap

            jj = 0
            do ii=1,length
            if(daily_OCO2(ii)%reject.eq.0) then
                jj=jj+1
                if(mod(jj,leap).ne.0) then
                daily_OCO2(ii)%reject = 1
                endif
            endif
            enddo

        endif

        totalnum = length - sum(daily_OCO2(:)%reject)
        write(*,*) "num in one day after thinning:    ",totalnum

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine OCOGST_QualityCtl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine OCOGST_Operator( obsdate, filename_in, length, Win_obs, DA_obs)
    !==========================================================================
    !!
    !!  OCO2 and gosat obs operator
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(8),intent(in)                 ::  obsdate
        character(*),intent(in)                 ::  filename_in
        !   filename_in pwd of TTRun
        integer,intent(in)                      ::  length
        type(Satellite_obs_bmc),intent(in)      ::  Win_obs(length)
        type(DA_obs_bmc),intent(out)            ::  DA_obs(length)

        !   loop
        integer             ::  ii,     jj,     kk

        !   inner var
        character(120)      ::  filename
        real(8)             ::  CO2_inst( nx, ny, nz, ntimes)
        real(8)             ::  PRES_inst( nx, ny, nz, ntimes)

        !------------------------------------------------------------------

        !   assign
        write(*,*) "CHECK CAREFULLY!"
        DA_obs(:)%lon       = Win_obs(:)%lon * 1.0d0
        DA_obs(:)%lat       = Win_obs(:)%lat * 1.0d0
        DA_obs(:)%alt       = -999.0d0
        DA_obs(:)%obs       = Win_obs(:)%xco2 * 1.0d0
        DA_obs(:)%obs_err   = Win_obs(:)%xco2_err * 1.0d0

        !   read in xb co2, and pres
        filename = trim(filename_in)//"/Backrun/DigInst."//obsdate//".bpch"
        call Extract_InstCO2( filename, tmp_CO2_r4, tmp_Prue_r4)
        CO2_inst(:,:,:,:)   = tmp_CO2_r4(:,:,:,:) * 1.0d0 * ppm
        PRES_inst(:,:,:,:)  = tmp_Prue_r4(:,:,:,:) * 1.0d0

        call OCOGST_Optcore( CO2_inst, PRES_inst, length, Win_obs, DA_obs(:)%obs_Xb)

        !   read in AA co2, and pres
        do ii=1,NRENS
        filename = trim(filename_in)//"/"//NumCha(ii)//"/DigInst."//obsdate//".bpch"
        call Extract_InstCO2( filename, tmp_CO2_r4, tmp_Prue_r4)
        CO2_inst(:,:,:,:)   = tmp_CO2_r4(:,:,:,:) * 1.0d0 * ppm
        PRES_inst(:,:,:,:)  = tmp_Prue_r4(:,:,:,:) * 1.0d0

        call OCOGST_Optcore( CO2_inst, PRES_inst, length, Win_obs, DA_obs(:)%obs_AA(ii))
        enddo

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine OCOGST_Operator

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine OCOGST_Optcore(  CO2_inst, PRES_inst, length, Win_obs, xco2_optd)
    !==========================================================================
    !!
    !!  OCO2 and gosat obs operator core from CO2 and pres to Simulated XCO2(:).
    !!
    !!  notice:
    !!      obs latitude range : -80.0 ~ 80.0
    !!      find star point is one good way to locate the obs point.
    !!
    !!      all at center:
    !!         0.0   2.5   5.0   7.5
    !!          ._____._____._____.  6.0
    !!          |     |     |     |
    !!          |     |     |     |
    !!          ._____._____._____.  4.0
    !!          |     |.....|     |
    !!          |     |.....|     |
    !!          ._____._____._____.  2.0
    !!          |     |     |     |
    !!          |     |     |     |
    !!          ._____._____._____.  0.0
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        real(8),intent(in)                  ::  CO2_inst( nx, ny, nz, ntimes)
        real(8),intent(in)                  ::  PRES_inst( nx, ny, nz, ntimes)
        integer,intent(in)                  ::  length
        type(Satellite_obs_bmc),intent(in)  ::  Win_obs(length)
        real(8),intent(out)                 ::  xco2_optd(length)

        !   loop
        integer                 ::  ii, jj, kk, ll, tt
        
        !   start point
        integer                 ::  stpt_lon,   stpt_lat

        !   value
        real(8)                 ::  CO2_col ( nbox, nbox, nz)
        real(8)                 ::  PRES_col( nbox, nbox, nz)
        !   dim_axis
        !   Point_lon rather than Point_lon
        !   plane_4lagr, must use vector.
        real(8)                 ::  Box_lon(nbox), Point_lon
        real(8)                 ::  Box_lat(nbox), Point_lat
        real(8)                 ::  Box_wght( nbox, nbox)

        !   results
        real(8)                 ::  CO2_optd (nz)
        real(8)                 ::  PRES_optd(nz)

        real(8)                 ::  CO2_prof_optd(OCOGST_level)
        real(8)                 ::  xco2_tmp

        !------------------------------------------------------------------

        do ll=1,length

            !   locate start point (checked)
            call locate_startpoint( Win_obs(ll)%lon, Win_obs(ll)%lat, stpt_lon, stpt_lat )

            !   locate value (checked)
            tt = int(Win_obs(ll)%hour/3) + 1
            call locate_Boxvalue( stpt_lon, stpt_lat, CO2_inst(:,:,:,tt), CO2_col)
            call locate_Boxvalue( stpt_lon, stpt_lat, PRES_inst(:,:,:,tt), PRES_col)

            !   locate axis (checked)
            call locate_Boxaxis( Win_obs(ll)%lon, Win_obs(ll)%lat, &
                &Box_lon, Box_lat, Point_lon, Point_lat )

            !   interpolation : locate weight (checked)
            call locate_weight(Box_lon, Box_lat, Point_lon, Point_lat, Box_wght)
            do ii=1,nz
                CO2_optd(ii)    = sum(Box_wght(:,:) * CO2_col(:,:,ii))
                PRES_optd(ii)   = sum(Box_wght(:,:) * PRES_col(:,:,ii))
            enddo

            !   vertical interpolation (checked)
            call line_4lagr( nz, CO2_optd(nz:1:-1), PRES_optd(nz:1:-1), &
                    &OCOGST_level, Win_obs(ll)%pres_lev(:)*1.0d0, CO2_prof_optd(:))

            !   compute xco2_optd (checked)
            xco2_optd(ll) = 0.0d0
            do ii=1,OCOGST_level
                xco2_tmp = CO2_prof_optd(ii) - Win_obs(ll)%co2_prof(ii)
                xco2_tmp = xco2_tmp * Win_obs(ll)%xco2_avk(ii) * Win_obs(ll)%pres_wgh(ii)
                xco2_optd(ll) = xco2_optd(ll) + xco2_tmp
            enddo
            xco2_optd(ll) = xco2_optd(ll) + Win_obs(ll)%xco2_apriori

        enddo

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine OCOGST_Optcore

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine locate_weight(Box_lon, Box_lat, Point_lon, Point_lat, Box_wght)
    !==========================================================================
    !!
    !!  locate Box val weight
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        real(8),intent(in)      ::  Box_lon(nbox)
        real(8),intent(in)      ::  Box_lat(nbox)
        real(8),intent(in)      ::  Point_lon
        real(8),intent(in)      ::  Point_lat
        real(8),intent(out)     ::  Box_wght(nbox,nbox)

        integer                 ::  ii, jj, kk

        integer                 ::  MM
        integer                 ::  NN

        !------------------------------------------------------------------
        !   
        !   lon1        lon2
        !   .-----------.
        !   MM          1-MM
        !
        !   lat2 . 1-NN
        !        |
        !        |
        !        |
        !   lat1 . NN

        MM =  ( Box_lon(2) - Point_lon ) / ( Box_lon(2) - Box_lon(1) )
        NN =  ( Box_lat(2) - Point_lat ) / ( Box_lat(2) - Box_lat(1) )

        Box_wght(1,1) = MM * NN
        Box_wght(1,2) = ( 1 - NN ) * MM
        Box_wght(2,1) = ( 1 - MM ) * NN
        Box_wght(2,2) = ( 1 - MM ) * ( 1 - NN )

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine locate_weight

    subroutine locate_Boxaxis( lon, lat, Box_lon, Box_lat, Point_lon, Point_lat )
    !==========================================================================
    !!
    !!  locate Box axis value(Box_lon, Box_lat) from Point_lon, Point_lat
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        real,intent(in)         ::  lon
        real,intent(in)         ::  lat
        real(8),intent(out)     ::  Box_lon(nbox)
        real(8),intent(out)     ::  Box_lat(nbox)
        real(8),intent(out)     ::  Point_lon
        real(8),intent(out)     ::  Point_lat

        integer                 ::  ii, jj, kk

        !------------------------------------------------------------------

        do ii=1,nbox
            Box_lon(ii) = (ii - 1) * 2.5d0
            Box_lat(ii) = (ii - 1) * 2.0d0
        enddo
        Point_lon = abs( floor( lon / 2.5d0 ) - lon / 2.5d0 ) * 2.5d0
        Point_lat = abs( floor( lat / 2.0d0 ) - lat / 2.0d0 ) * 2.0d0

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine locate_Boxaxis

    subroutine locate_Boxvalue(  stpt_lon, stpt_lat, Var_inst, Var_col)
    !==========================================================================
    !!
    !!  locate Box column value(Var_col) from DigInst value(Var_inst)
    !!  with start point (stpt_lon,stpt_lat)
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        integer,intent(in)      ::  stpt_lon
        integer,intent(in)      ::  stpt_lat
        real(8),intent(in)      ::  Var_inst( nx, ny, nz)
        real(8),intent(out)     ::  Var_col ( nbox, nbox, nz)

        integer                 ::  ii, jj, kk

        !------------------------------------------------------------------
        !   how to deal with abnormal stpt_lon
        !   nx-3   nx-2   nx-1     nx      1      2      3      4
        !      .      .      .      .      .      .      .      .
        !      .      .      .      .
        !
        !             *      .      .      .
        !                    *      .      .      .
        !                           *      .      .      .
        !
        !                                  .      .      .      .
        !

        if( stpt_lon.lt.nx .and. stpt_lon.ge.1 ) then
            Var_col(:,:,:)  = Var_inst( stpt_lon:stpt_lon+1, stpt_lat:stpt_lat+1, :)
        elseif( stpt_lon .eq. nx ) then
            Var_col(1,:,:)  = Var_inst( nx, stpt_lat:stpt_lat+1, :)
            Var_col(2,:,:)  = Var_inst( 1,  stpt_lat:stpt_lat+1, :)
        else
            write(*,*) stpt_lon
            stop "check locate_Boxvalue"
        endif

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine locate_Boxvalue

    subroutine locate_startpoint( lon, lat, stpt_lon, stpt_lat)
    !==========================================================================
    !!
    !!  locate start point of lon and lat
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        real,intent(in)         ::  lon 
        real,intent(in)         ::  lat
        integer,intent(out)     ::  stpt_lon
        integer,intent(out)     ::  stpt_lat

        integer                 ::  ii, jj, kk

        !------------------------------------------------------------------

        !   locate stpt_lon
        stpt_lon = 0
        if(     lon*1.0d0 .ge. MLonC(nx) .and. lon*1.0d0 .le. 180.0d0     ) then
        stpt_lon = nx
        else
        loop1: do ii=1,nx-1
            if( lon*1.0d0 .ge. MLonC(ii) .and. lon*1.0d0 .lt. MLonC(ii+1) ) then
            stpt_lon = ii
            exit loop1
            endif
        enddo loop1
        endif

        !   locate stpt_lat
        stpt_lat = 0
        loop2: do jj=1,ny-1
            if( lat*1.0d0 .gt. MLatC(jj) .and. lat*1.0d0 .le. MLatC(jj+1) ) then
            stpt_lat = jj
            exit loop2
            endif
        enddo loop2

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine locate_startpoint

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine handle_err( Status_err, filename_in, AName_err)
    !==========================================================================
    !!
    !!  handle error
    !!
    !==========================================================================
        implicit none
        !------------------------------------------------------------------

        integer,intent(in)              ::  Status_err
        character(*),intent(in)         ::  filename_in
        character(*),intent(in)         ::  AName_err

        if (Status_err .ne. NF_NOERR) then
            print *, trim(filename_in)
            print *, trim(AName_err)
            print *, nf_strerror(Status_err)
            stop 'stopped'
        end if

        return
        !------------------------------------------------------------------

    !==========================================================================
    end subroutine handle_err

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==================================================================================
end MODULE operator_ocogst
