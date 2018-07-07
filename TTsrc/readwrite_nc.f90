MODULE ReadWrite_nc
!!=================================================================================
!!
!!  all the operations are in R4!
!!
!!    Read and Write netcdf format file of
!!      CoeRe(nx, ny, nCoeRe, 1)
!!      AREA(nx, ny)
!!      CO2(nx, ny, nz, 1)
!!      Flux(nx, ny, nCoeRe, 1)
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2017.12.27                   Hanrui                 NV2.0 original code
!!     2018.03.28                   Hanrui                 NV2.0 Final version
!!
!!=================================================================================

    !   almost all variables in Common are needed
    use Common
    use netcdf
    implicit none
    include 'netcdf.inc'
    !--------------------------------------------------------------------------
    !   subroutine titles
    !   CoeRe
    PUBLIC                  ::  Write_CoeRe_nc
    PUBLIC                  ::  Read_CoeRe_nc
    !   CoeRe
    PUBLIC                  ::  Write_Flux_nc
    PUBLIC                  ::  Read_Flux_nc
    !   grid area
    PRIVATE                 ::  Write_AREA_nc
    PUBLIC                  ::  Read_AREA_nc
    !   CO2
    PUBLIC                  ::  Write_CO2_nc
    PUBLIC                  ::  Read_CO2_nc
    !read from r4 to r8
    real(4),PUBLIC          ::  tmp_CoeRe_r4(nx,ny,nCoeRe,1)
    real(4),PUBLIC          ::  tmp_Flux_r4(nx,ny,nCoeRe,1)
    real(4),PUBLIC          ::  tmp_CO2nc_r4(nx,ny,nz,1)
    real(8),PUBLIC          ::  AREAgird(nx,ny)

    !--------------------------------------------------------------------------
    PRIVATE

        !------------------------------------------------------------------
        !inner var:
        !   Attributes char
        character(len=256)      ::  AChar
        character(len=256)      ::  AName
        !   file id and Status
        integer                 ::  Status
        integer                 ::  ncid
        !   coordinate dimension id and variable id
        integer                 ::  x_dimid,    x_varid
        integer                 ::  y_dimid,    y_varid
        integer                 ::  z_dimid,    z_varid
        integer                 ::  c_dimid,    c_varid
        integer                 ::  t_dimid,    t_varid
        !   variable id (all of them are 4 dimensions)
        integer,dimension(4)    ::  dimid
        integer                 ::  varid
        !   gird area
        integer,dimension(2)    ::  ga_dimid
        integer                 ::  ga_varid
        !   time char
        character(10)           ::  TimeCha

        !------------------------------------------------------------------

    !--------------------------------------------------------------------------
    Contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Read_CoeRe_nc( infile, CoeRe)
    !==========================================================================
    !!
    !!  Read in CoeRe(nx,ny,nCoeRe,1) from netcdf file
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !outputs:
        !   CoeRe
        real,intent(out)        ::  CoeRe(nx,ny,nCoeRe,1)

        !------------------------------------------------------------------

        !OPen NETCDF file
        Status = nf_open( trim(infile), nf_clobber, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) 'OPEN ERROR CHECK_FILE: ', trim(infile)
        stop
        end if

        !inquire varid
        AName  = 'CoeRe'
        Status = nf_inq_varid(ncid, trim(AName), varid)
        if (Status /= NF_NOERR) then
        write(*,*) 'INQUIRE ERROR CHECK_FILE: ', trim(infile)
        write(*,*) 'INQUIRE ERROR CHECK_VAR: ',  trim(AName)
        stop
        end if

        !read var
        Status = nf_get_var_real(ncid, varid, CoeRe)

        !close
        Status = nf_close( ncid)

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Read_CoeRe_nc

    subroutine Write_CoeRe_nc( infile, DateCha, CoeRe)
    !==========================================================================
    !!
    !!  Write Out CoeRe(nx,ny,nCoeRe,1) into netcdf file
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !   CoeRe time
        character(8),intent(in) ::  DateCha
        !   CoeRe
        real,intent(in)         ::  CoeRe(nx,ny,nCoeRe,1)
        integer                 ::  ll

        !Set file time
        !------------------------------------------------------------------
        TimeCha(1:4)    = DateCha(1:4)
        TimeCha(5:5)    = "-"
        TimeCha(6:7)    = DateCha(5:6)
        TimeCha(8:8)    = "-"
        TimeCha(9:10)   = DateCha(7:8)

        !Create netCDF file
        !------------------------------------------------------------------
        !Create a new nc file.(Overwrite the existing one)
        write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        write(*,*) "Creating ",trim(infile)," ..."
        write(*,*) ""

        Status = nf_create( trim(infile), nf_clobber, ncid)

        !***

        !define the dimensions
        !------------------------------------------------------------------
        !Define dimensions
        !    In  : ncid, dim_name, dim_lenth
        !    Out : dim_id
        write(*,*) "Defining dimensions ..."

        Status = nf_def_dim( ncid, "lon",    nx,     x_dimid)
        Status = nf_def_dim( ncid, "lat",    ny,     y_dimid)
        Status = nf_def_dim( ncid, "clev",   nCoeRe, c_dimid)
        Status = nf_def_dim( ncid, "time",   1,      t_dimid)

        !***

        !define coordinate variables
        !------------------------------------------------------------------
        !Define coordinate variables of dimensions
        !    data_type needed here.
        !    In  :ncid, dim_name, data_type, 1(dimension), dim_id
        !    Out :var id
        write(*,*) "Defining coordinate variables of dimensions ..."

        Status = nf_def_var( ncid, "lon",   nf_double, 1, x_dimid, x_varid)
        Status = nf_def_var( ncid, "lat",   nf_double, 1, y_dimid, y_varid)
        Status = nf_def_var( ncid, "clev",  nf_double, 1, c_dimid, c_varid)
        Status = nf_def_var( ncid, "time",  nf_float,  1, t_dimid, t_varid)

        !define coordinate variables' attributes
        !------------------------------------------------------------------
        !Defining coordinate variables' attributes
        !   coordinate variables always have long_name and units only
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining coordinate variables' attributes ..."

        write(*,*) "    long_name ..."
        !   long_name
        AName  = "long_name"
        AChar  = "Longitude"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Latitude"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Flux species"
        Status = nf_put_att_text( ncid, c_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Time"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        !   units
        AName  = "units"
        AChar  = "degrees_east"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "degrees_north"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "unitless"
        Status = nf_put_att_text( ncid, c_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "days since "//TimeCha//" 00:00:00 GMT"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        !***

        !define variables
        !------------------------------------------------------------------
        !Define variables   
        write(*,*) ""
        write(*,*) "Defining variable CoeRe ... "

        write(*,*) "    Give dimension ..."
        dimid = (/ x_dimid, y_dimid, c_dimid, t_dimid /)

        write(*,*) "    Defining variable ..."
        Status = nf_def_var( ncid, "CoeRe", nf_float, 4, dimid, varid)

        !define variables' attributes
        !------------------------------------------------------------------
        !Defining variables' attributes
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining variable attributes ... "

        write(*,*) "    long_name ..."
        AName  = "long_name"
        AChar  = "GEOS-Chem 10 Species Flux Scale factor"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        AName  = "units"
        AChar  = "unitless scale factor after range ctl."
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    averaging_method ..."
        AName  = "averaging_method"
        AChar  = "Daily average"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    FillValue ..."
        AName  = "_FillValue"
        Status = nf_put_att_real( ncid, varid, trim(AName), nf_float, 1, -1.0E31)

        write(*,*) "    clev_name ..."
        do ll = 1, nCoeRe
        AName  = "clev_name_"//NumCha(ll)
        AChar  = FluxSpi(ll)
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))
        enddo

        !***

        !end defination
        !------------------------------------------------------------------
        write(*,*) ""
        write(*,*) "End define"
        Status = nf_enddef( ncid)

        !***

        !write data
        !------------------------------------------------------------------
        write(*,*) "Write data"
        Status = nf_put_var_double( ncid, x_varid, MLonC)
        Status = nf_put_var_double( ncid, y_varid, MLatC)
        Status = nf_put_var_double( ncid, c_varid, MCoeReC)
        Status = nf_put_var_real( ncid, t_varid, 1.0d0)
        Status = nf_put_var_real( ncid, varid, CoeRe)

        !***

        !Close
        !------------------------------------------------------------------
        Status = nf_close( ncid)

        !***
        
        return
        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Write_CoeRe_nc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Read_Flux_nc( infile, Flux)
    !==========================================================================
    !!
    !!  Read in Flux(nx,ny,nCoeRe,1) from netcdf file
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !outputs:
        !   Flux
        real,intent(out)        ::  Flux(nx,ny,nCoeRe,1)

        !------------------------------------------------------------------

        !OPen NETCDF file
        Status = nf_open( trim(infile), nf_clobber, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) 'OPEN ERROR CHECK_FILE: ', trim(infile)
        stop
        end if

        !inquire varid
        AName  = 'Flux'
        Status = nf_inq_varid(ncid, trim(AName), varid)
        if (Status /= NF_NOERR) then
        write(*,*) 'INQUIRE ERROR CHECK_FILE: ', trim(infile)
        write(*,*) 'INQUIRE ERROR CHECK_VAR: ',  trim(AName)
        stop
        end if

        !read var
        Status = nf_get_var_real(ncid, varid, Flux)

        !close
        Status = nf_close( ncid)

        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Read_Flux_nc


    subroutine Write_Flux_nc( infile, DateCha, Flux)
    !==========================================================================
    !!
    !!  Write Out Flux(nx,ny,nCoeRe,1) into netcdf file
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !   Flux time
        character(8),intent(in) ::  DateCha
        !   Flux
        real,intent(in)         ::  Flux(nx,ny,nCoeRe,1)
        integer                 ::  ll

        !Set file time
        !------------------------------------------------------------------
        TimeCha(1:4)    = DateCha(1:4)
        TimeCha(5:5)    = "-"
        TimeCha(6:7)    = DateCha(5:6)
        TimeCha(8:8)    = "-"
        TimeCha(9:10)   = DateCha(7:8)

        !Create netCDF file
        !------------------------------------------------------------------
        !Create a new nc file.(Overwrite the existing one)
        write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        write(*,*) "Creating ",trim(infile)," ..."
        write(*,*) ""

        Status = nf_create( trim(infile), nf_clobber, ncid)

        !***

        !define the dimensions
        !------------------------------------------------------------------
        !Define dimensions
        !    In  : ncid, dim_name, dim_lenth
        !    Out : dim_id
        write(*,*) "Defining dimensions ..."

        Status = nf_def_dim( ncid, "lon",    nx,     x_dimid)
        Status = nf_def_dim( ncid, "lat",    ny,     y_dimid)
        Status = nf_def_dim( ncid, "clev",   nCoeRe, c_dimid)
        Status = nf_def_dim( ncid, "time",   1,      t_dimid)

        !***

        !define coordinate variables
        !------------------------------------------------------------------
        !Define coordinate variables of dimensions
        !    data_type needed here.
        !    In  :ncid, dim_name, data_type, 1(dimension), dim_id
        !    Out :var id
        write(*,*) "Defining coordinate variables of dimensions ..."

        Status = nf_def_var( ncid, "lon",   nf_double, 1, x_dimid, x_varid)
        Status = nf_def_var( ncid, "lat",   nf_double, 1, y_dimid, y_varid)
        Status = nf_def_var( ncid, "clev",  nf_double, 1, c_dimid, c_varid)
        Status = nf_def_var( ncid, "time",  nf_float,  1, t_dimid, t_varid)

        !define coordinate variables' attributes
        !------------------------------------------------------------------
        !Defining coordinate variables' attributes
        !   coordinate variables always have long_name and units only
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining coordinate variables' attributes ..."

        write(*,*) "    long_name ..."
        !   long_name
        AName  = "long_name"
        AChar  = "Longitude"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Latitude"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Flux species"
        Status = nf_put_att_text( ncid, c_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Time"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        !   units
        AName  = "units"
        AChar  = "degrees_east"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "degrees_north"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "unitless"
        Status = nf_put_att_text( ncid, c_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "days since "//TimeCha//" 00:00:00 GMT"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        !***

        !define variables
        !------------------------------------------------------------------
        !Define variables   
        write(*,*) ""
        write(*,*) "Defining variable Flux ... "

        write(*,*) "    Give dimension ..."
        dimid = (/ x_dimid, y_dimid, c_dimid, t_dimid /)

        write(*,*) "    Defining variable ..."
        Status = nf_def_var( ncid, "Flux", nf_float, 4, dimid, varid)

        !define variables' attributes
        !------------------------------------------------------------------
        !Defining variables' attributes
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining variable attributes ... "

        write(*,*) "    long_name ..."
        AName  = "long_name"
        AChar  = "GEOS-Chem 10 Species Flux"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        AName  = "units"
        AChar  = "g C/m2/day"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    averaging_method ..."
        AName  = "averaging_method"
        AChar  = "Daily average"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    FillValue ..."
        AName  = "_FillValue"
        Status = nf_put_att_real( ncid, varid, trim(AName), nf_float, 1, -1.0E31)

        write(*,*) "    clev_name ..."
        do ll = 1, nCoeRe
        AName  = "clev_name_"//NumCha(ll)
        AChar  = FluxSpi(ll)
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))
        enddo

        !***

        !end defination
        !------------------------------------------------------------------
        write(*,*) ""
        write(*,*) "End define"
        Status = nf_enddef( ncid)

        !***

        !write data
        !------------------------------------------------------------------
        write(*,*) "Write data"
        Status = nf_put_var_double( ncid, x_varid, MLonC)
        Status = nf_put_var_double( ncid, y_varid, MLatC)
        Status = nf_put_var_double( ncid, c_varid, MCoeReC)
        Status = nf_put_var_real( ncid, t_varid, 1.0d0)
        Status = nf_put_var_real( ncid, varid, Flux)

        !***

        !Close
        !------------------------------------------------------------------
        Status = nf_close( ncid)

        !***
        
        return
        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Write_Flux_nc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Read_AREA_nc( infile, AREA)
    !==========================================================================
    !!
    !!  Read in AREA(nx,ny) from netcdf file
    !!  
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !outputs:
        !   AREA
        real(8),intent(out)     ::  AREA(nx,ny)

        !------------------------------------------------------------------

        !OPen NETCDF file
        Status = nf_open( trim(infile), nf_clobber, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) 'OPEN ERROR CHECK_FILE: ', trim(infile)
        stop
        end if

        !inquire varid
        AName  = "DXYP"
        Status = nf_inq_varid(ncid, trim(AName), varid)
        if (Status /= NF_NOERR) then
        write(*,*) 'INQUIRE ERROR CHECK_FILE: ', trim(infile)
        write(*,*) 'INQUIRE ERROR CHECK_VAR: ',  trim(AName)
        stop
        end if

        !read var
        Status = nf_get_var_double(ncid, varid, AREA)

        !close
        Status = nf_close( ncid)
        
        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Read_AREA_nc

    subroutine Write_AREA_nc( infile, AREA)
    !==========================================================================
    !!
    !!  Write Out AREA(nx,ny) into netcdf file (use only once)
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !   AREA
        real(8),intent(in)      ::  AREA(nx,ny)

        !------------------------------------------------------------------

        !Create netCDF file
        !------------------------------------------------------------------
        !Create a new nc file.(Overwrite the existing one)
        write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        write(*,*) "Creating ",trim(infile)," ..."
        write(*,*) ""

        Status = nf_create( trim(infile), nf_clobber, ncid)

        !***

        !define the dimensions
        !------------------------------------------------------------------
        !Define dimensions
        !    In  : ncid, dim_name, dim_lenth
        !    Out : dim_id
        write(*,*) "Defining dimensions ..."

        Status = nf_def_dim( ncid, "lon",    nx,     x_dimid)
        Status = nf_def_dim( ncid, "lat",    ny,     y_dimid)

        !***

        !define coordinate variables
        !------------------------------------------------------------------
        !Define coordinate variables of dimensions
        !    data_type needed here.
        !    In  :ncid, dim_name, data_type, 1(dimension), dim_id
        !    Out :var id
        write(*,*) "Defining coordinate variables of dimensions ..."

        Status = nf_def_var( ncid, "lon",   nf_double, 1, x_dimid, x_varid)
        Status = nf_def_var( ncid, "lat",   nf_double, 1, y_dimid, y_varid)

        !define coordinate variables' attributes
        !------------------------------------------------------------------
        !Defining coordinate variables' attributes
        !   coordinate variables always have long_name and units only
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining coordinate variables' attributes ..."

        write(*,*) "    long_name ..."
        !   long_name
        AName  = "long_name"
        AChar  = "Longitude"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Latitude"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        !   units
        AName  = "units"
        AChar  = "degrees_east"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "degrees_north"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))

        !***

        !define variables
        !------------------------------------------------------------------
        !Define variables   
        write(*,*) ""
        write(*,*) "Defining variable CoeRe ... "

        write(*,*) "    Give dimension ..."
        ga_dimid = (/ x_dimid, y_dimid /)

        write(*,*) "    Defining variable ..."
        Status = nf_def_var( ncid, "DXYP", nf_double, 2, ga_dimid, ga_varid)

        !define variables' attributes
        !------------------------------------------------------------------
        !Defining variables' attributes
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining variable attributes ... "

        write(*,*) "    long_name ..."
        AName  = "long_name"
        AChar  = "Grid box area"
        Status = nf_put_att_text( ncid, ga_varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        AName  = "units"
        AChar  = "m2"
        Status = nf_put_att_text( ncid, ga_varid, trim(AName), len_trim(AChar), trim(AChar))

        !***

        !end defination
        !------------------------------------------------------------------
        write(*,*) ""
        write(*,*) "End define"
        Status = nf_enddef( ncid)

        !***

        !write data
        !------------------------------------------------------------------
        write(*,*) "Write data"
        Status = nf_put_var_double( ncid, x_varid, MLonC)
        Status = nf_put_var_double( ncid, y_varid, MLatC)
        Status = nf_put_var_double( ncid, ga_varid, AREA)

        !***

        !Close
        !------------------------------------------------------------------
        Status = nf_close( ncid)

        !***
        
        return
        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Write_AREA_nc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Read_CO2_nc( infile, CO2)
    !==========================================================================
    !!
    !!  Read in CO2(nx,ny,nz,1) from netcdf file
    !!  
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        !outputs:
        !   CoeRe
        real,intent(out)        ::  CO2(nx,ny,nz,1)

        !------------------------------------------------------------------

        !OPen NETCDF file
        Status = nf_open( trim(infile), nf_clobber, ncid)
        if (Status /= NF_NOERR) then
        write(*,*) 'OPEN ERROR CHECK_FILE: ', trim(infile)
        stop
        end if

        !inquire varid
        AName  = "SPC_CO2"
        Status = nf_inq_varid(ncid, trim(AName), varid)
        if (Status /= NF_NOERR) then
        write(*,*) 'INQUIRE ERROR CHECK_FILE: ', trim(infile)
        write(*,*) 'INQUIRE ERROR CHECK_VAR: ',  trim(AName)
        stop
        end if

        !read var
        Status = nf_get_var_real(ncid, varid, CO2)

        !close
        Status = nf_close( ncid)
        
        return

        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Read_CO2_nc

    subroutine Write_CO2_nc( infile, DateCha, CO2)
    !==========================================================================
    !!
    !!  Write Out CO2(nx,ny,nz,1) into netcdf file
    !!
    !==========================================================================

        implicit none

        !------------------------------------------------------------------
        !inputs:
        !   file name
        character(*),intent(in) ::  infile
        character(8),intent(in) ::  DateCha
        !   CO2
        real,intent(in)         ::  CO2(nx,ny,nz,1)

        !Set file time
        !------------------------------------------------------------------
        TimeCha(1:4)    = DateCha(1:4)
        TimeCha(5:5)    = "-"
        TimeCha(6:7)    = DateCha(5:6)
        TimeCha(8:8)    = "-"
        TimeCha(9:10)   = DateCha(7:8)

        !Create netCDF file
        !------------------------------------------------------------------
        !Create a new nc file.(Overwrite the existing one)
        write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
        write(*,*) "Creating ",trim(infile)," ..."
        write(*,*) ""

        Status = nf_create( trim(infile), nf_clobber, ncid)

        !***

        !define the dimensions
        !------------------------------------------------------------------
        !Define dimensions
        !    In  : ncid, dim_name, dim_lenth
        !    Out : dim_id
        write(*,*) "Defining dimensions ..."

        Status = nf_def_dim( ncid, "lon",    nx,     x_dimid)
        Status = nf_def_dim( ncid, "lat",    ny,     y_dimid)
        Status = nf_def_dim( ncid, "lev",    nz,     z_dimid)
        Status = nf_def_dim( ncid, "time",   1,      t_dimid)

        !***

        !define coordinate variables
        !------------------------------------------------------------------
        !Define coordinate variables of dimensions
        !    data_type needed here.
        !    In  :ncid, dim_name, data_type, 1(dimension), dim_id
        !    Out :var id
        write(*,*) "Defining coordinate variables of dimensions ..."

        Status = nf_def_var( ncid, "lon",   nf_double, 1, x_dimid, x_varid)
        Status = nf_def_var( ncid, "lat",   nf_double, 1, y_dimid, y_varid)
        Status = nf_def_var( ncid, "lev",   nf_double, 1, z_dimid, z_varid)
        Status = nf_def_var( ncid, "time",  nf_float,  1, t_dimid, t_varid)

        !define coordinate variables' attributes
        !------------------------------------------------------------------
        !Defining coordinate variables' attributes
        !   coordinate variables always have long_name and units only
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining coordinate variables' attributes ..."

        write(*,*) "    long_name ..."
        !   long_name
        AName  = "long_name"
        AChar  = "Longitude"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Latitude"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "GEOS-Chem level"
        Status = nf_put_att_text( ncid, z_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "Time"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        !   units
        AName  = "units"
        AChar  = "degrees_east"
        Status = nf_put_att_text( ncid, x_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "degrees_north"
        Status = nf_put_att_text( ncid, y_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "unitless"
        Status = nf_put_att_text( ncid, z_varid, trim(AName), len_trim(AChar), trim(AChar))
        AChar  = "hours since "//TimeCha//" 00:00:00 GMT"
        Status = nf_put_att_text( ncid, t_varid, trim(AName), len_trim(AChar), trim(AChar))

        !***

        !define variables
        !------------------------------------------------------------------
        !Define variables   
        write(*,*) ""
        write(*,*) "Defining variable SPC_CO2 ... "

        write(*,*) "    Give dimension ..."
        dimid = (/ x_dimid, y_dimid, z_dimid, t_dimid /)

        write(*,*) "    Defining variable ..."
        Status = nf_def_var( ncid, "SPC_CO2", nf_float, 4, dimid, varid)

        !define variables' attributes
        !------------------------------------------------------------------
        !Defining variables' attributes
        !   data_type needed here.
        !   In  :ncid, var id, attribute_name, dimensions, variable
        !   Out :
        write(*,*) "Defining variable attributes ... "

        write(*,*) "    long_name ..."
        AName  = "long_name"
        AChar  = "SPC_CO2"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    units ..."
        AName  = "units"
        AChar  = "mol mol-1"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    averaging_method ..."
        AName  = "averaging_method"
        AChar  = "Instantaneous"
        Status = nf_put_att_text( ncid, varid, trim(AName), len_trim(AChar), trim(AChar))

        write(*,*) "    FillValue ..."
        AName  = "_FillValue"
        Status = nf_put_att_real( ncid, varid, trim(AName), nf_float, 1, -1.0E31)
        
        !***

        !end defination
        !------------------------------------------------------------------
        write(*,*) ""
        write(*,*) "End define"
        Status = nf_enddef( ncid)

        !***

        !write data
        !------------------------------------------------------------------
        write(*,*) "Write data"
        Status = nf_put_var_double( ncid, x_varid, MLonC)
        Status = nf_put_var_double( ncid, y_varid, MLatC)
        Status = nf_put_var_double( ncid, z_varid, MLevC)
        Status = nf_put_var_real( ncid, varid, CO2)
        Status = nf_put_var_real( ncid, t_varid, 1.0d0)

        !***

        !Close
        !------------------------------------------------------------------
        Status = nf_close( ncid)

        !***
        
        return
        !------------------------------------------------------------------

    !==========================================================================
    end subroutine Write_CO2_nc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !==========================================================================
    
!!=================================================================================
END MODULE ReadWrite_nc