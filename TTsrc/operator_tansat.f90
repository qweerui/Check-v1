MODULE operator_tansat
!!=================================================================================
!!
!!  Description:
!!      OBS operator of tansat.
!!
!!  Warning:
!!      var(lon,lat,lev) in nc_file, into var(lev,lat,lon) in fortran.
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2018.01.02                   Hanrui                  NV2.0 original code
!!
!!=================================================================================

    use common
    use netcdf

    use InputOption,    only :  NRENS
    use InputOption,    only :  TANSAT_level
    use InputOption,    only :  maxval_sat
    use InputOption,    only :  TanSat_obs_bmc
    use InputOption,    only :  DA_obs_bmc

    use interpolation,  only :  plane_4lagr,line_4lagr
    use interpolation,  only :  plane_3lagr,line_3lagr
    use interpolation,  only :  plane_2lagr,line_2lagr

    use extract_bpch,   only :  Extract_InstCO2
    use extract_bpch,   only :  tmp_CO2_r4,tmp_Prue_r4
    
    use operator_ocogst,only :  locate_startpoint
    use operator_ocogst,only :  locate_Boxvalue
    use operator_ocogst,only :  locate_Boxaxis
    use operator_ocogst,only :  locate_weight
    use operator_ocogst,only :  handle_err


    implicit none
    include "netcdf.inc"
    !--------------------------------------------------------------------------
    !   find TANSAT or gosat file name and length
    PUBLIC                  ::  TANSAT_Filename_length
    !   read in TANSAT or gosat obs
    PUBLIC                  ::  TANSAT_read
    !   TANSAT obs quality ctl
    PUBLIC                  ::  TANSAT_QualityCtl
    !   TANSAT obs operator
    PUBLIC                  ::  TANSAT_Operator
    !   TANSAT obs operator from CO2 and pres to Simulated XCO2(:)
    PUBLIC                  ::  TANSAT_Optcore



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

    subroutine TANSAT_Filename_length( obsdate, filename_in, length)
    !==========================================================================
    !!
    !!  find TANSAT file name and length
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(8),intent(in)         ::  obsdate
        character(*),intent(inout)      ::  filename_in
        integer,intent(out)             ::  length


        character(120)                  ::  namechar
        integer                         ::  Status
        !   loop
        integer                         ::  ii,     jj,     kk

        !------------------------------------------------------------------
        filename_in = trim(filename_in)//"/Obs/TanSat/04/TanSat_L2SCI_"
        filename_in = trim(filename_in)//obsdate(7:8)//".nc"


        !OPen NETCDF file
        Status = nf_open( trim(filename_in), NF_NOWRITE, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) "OPEN ERROR CHECK_FILE: ", trim(filename_in)
        stop
        end if

        !inquire varid dimension
        AName  = "c"
        Status = nf_inq_dimid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)

        !read var
        Status = nf_inq_dimlen( ncid, varid, length)
        call handle_err( Status, filename_in, AName)

        !close
        Status = nf_close( ncid)

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine TANSAT_Filename_length

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine TANSAT_read( filename_in, length, daily_TanSat)
    !==========================================================================
    !!
    !!  read in TANSAT obs
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(*),intent(in)                 ::  filename_in
        integer,intent(in)                      ::  length
        type(TanSat_obs_bmc),intent(out)        ::  daily_TanSat(length)

        !------------------------------------------------------------------
        !   
        real(8)             ::  double0D( 1, TANSAT_level )
        real(8)             ::  double1D( 1, length )
        real(8)             ::  double2D( TANSAT_level, length )
        real(8)             ::  double3D( TANSAT_level, TANSAT_level, length )
        !   loop
        integer             ::  ii,     jj,     kk

        !------------------------------------------------------------------

        !OPen NETCDF file
        write(*,*) trim(filename_in)
        Status = nf_open( trim(filename_in), NF_NOWRITE, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) "OPEN ERROR CHECK_FILE: ", trim(filename_in)
        stop
        end if

        !read year
        AName  = "year"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%year = int(double1D(1,ii))
        enddo

        !read month
        AName  = "month"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%month = int(double1D(1,ii))
        enddo

        !read day
        AName  = "day"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%day = int(double1D(1,ii))
        enddo

        !read hour
        AName  = "hour"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%hour = int(double1D(1,ii))
        enddo

        !read longitude
        AName  = "longitude"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%lon = real(double1D(1,ii))
        enddo

        !read latitude
        AName  = "latitude"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%lat = real(double1D(1,ii))
        enddo

        !read xco2_bc
        AName  = "xco2_bc"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%xco2 = real(double1D(1,ii))
        enddo
        
        !read xco2_error
        AName  = "xco2_error"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%xco2_err = real(double1D(1,ii))
        enddo

        !read xco2_apriori
        AName  = "xco2_apriori"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double1D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%xco2_apriori = real(double1D(1,ii))
        enddo

        !read xco2_column_averaging_kernel
        AName  = "xco2_column_averaging_kernel"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double2D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%xco2_avk(:) = real(double2D(:,ii))
        enddo

        !read co2_profile_apriori
        AName  = "co2_profile_apriori"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double2D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%co2_prof(:) = real(double2D(:,ii))
        enddo

        !read h presure weight
        AName  = "h"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double0D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%pres_wgh(:) = real(double0D(1,:))
        enddo

        !read pressure_levels
        AName  = "pressure_levels"
        Status = nf_inq_varid( ncid, trim(AName), varid)
        call handle_err( Status, filename_in, AName)
        Status = nf_get_var_double( ncid, varid, double2D)
        call handle_err( Status, filename_in, AName)
        do ii=1,length
        daily_TanSat(ii)%pres_lev(:) = real(double2D(:,ii))
        enddo

        ! write(*,*) daily_TanSat(1)
        ! write(*,*) daily_TanSat(length)
        ! stop

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine TANSAT_read

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine TANSAT_QualityCtl( length, daily_TanSat )
    !==========================================================================
    !!
    !!  OCO2 obs quality ctl
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        integer,intent(in)                      ::  length
        type(TanSat_obs_bmc),intent(inout)      ::  daily_TanSat(length)
        !   loop
        integer         ::  ii,     jj,     kk

        integer         ::  totalnum
        integer         ::  leap

        !------------------------------------------------------------------

        !   reject = 0 : accept
        !   reject = 1 : reject
        daily_TanSat(:)%reject = 0
        totalnum = length - sum(daily_TanSat(:)%reject)
        write(*,*) "Total num in one day:             ",totalnum

        !   data quality ctl
        do ii=1,length
        if(daily_TanSat(ii)%xco2_err.gt.3.0) daily_TanSat(ii)%reject = 1
        if(daily_TanSat(ii)%lat.gt.80.0 .or. daily_TanSat(ii)%lat.lt.-80.0) then
            daily_TanSat(ii)%reject = 1
        endif
        enddo

        totalnum = length - sum(daily_TanSat(:)%reject)
        write(*,*) "num in one day after quality ctl: ",totalnum

        
        !   data thinning
        if(totalnum.gt.maxval_sat) then

            leap = ceiling((totalnum*1.0)/(maxval_sat*1.0))
            !   ceiling : int upword
            write(*,*) "thinning leap use:                ",leap

            jj = 0
            do ii=1,length
            if(daily_TanSat(ii)%reject.eq.0) then
                jj=jj+1
                if(mod(jj,leap).ne.0) then
                daily_TanSat(ii)%reject = 1
                endif
            endif
            enddo

        endif

        totalnum = length - sum(daily_TanSat(:)%reject)
        write(*,*) "num in one day after thinning:    ",totalnum

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine TANSAT_QualityCtl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine TANSAT_Operator( obsdate, filename_in, length, Win_obs, DA_obs)
    !==========================================================================
    !!
    !!  TANSAT obs operator
    !!
    !==========================================================================

        implicit none
        !------------------------------------------------------------------
        character(8),intent(in)                 ::  obsdate
        character(*),intent(in)                 ::  filename_in
        !   filename_in pwd of TTRun
        integer,intent(in)                      ::  length
        type(TanSat_obs_bmc),intent(in)         ::  Win_obs(length)
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

        call TANSAT_Optcore( CO2_inst, PRES_inst, length, Win_obs, DA_obs(:)%obs_Xb)

        !   read in AA co2, and pres
        do ii=1,NRENS
        filename = trim(filename_in)//"/"//NumCha(ii)//"/DigInst."//obsdate//".bpch"
        call Extract_InstCO2( filename, tmp_CO2_r4, tmp_Prue_r4)
        CO2_inst(:,:,:,:)   = tmp_CO2_r4(:,:,:,:) * 1.0d0 * ppm
        PRES_inst(:,:,:,:)  = tmp_Prue_r4(:,:,:,:) * 1.0d0

        call TANSAT_Optcore( CO2_inst, PRES_inst, length, Win_obs, DA_obs(:)%obs_AA(ii))
        enddo

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine TANSAT_Operator

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine TANSAT_Optcore(  CO2_inst, PRES_inst, length, Win_obs, xco2_optd)
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
        type(TanSat_obs_bmc),intent(in)     ::  Win_obs(length)
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

        real(8)                 ::  CO2_prof_optd(TANSAT_level)
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
                    &TANSAT_level, Win_obs(ll)%pres_lev(:)*1.0d0, CO2_prof_optd(:))

            !   compute xco2_optd (checked)
            xco2_optd(ll) = 0.0d0
            do ii=1,TANSAT_level
                xco2_tmp = CO2_prof_optd(ii) - Win_obs(ll)%co2_prof(ii)
                xco2_tmp = xco2_tmp * Win_obs(ll)%xco2_avk(ii) * Win_obs(ll)%pres_wgh(ii)
                xco2_optd(ll) = xco2_optd(ll) + xco2_tmp
            enddo
            xco2_optd(ll) = xco2_optd(ll) + Win_obs(ll)%xco2_apriori

        enddo

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine TANSAT_Optcore

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!==================================================================================
end MODULE operator_tansat
