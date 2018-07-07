Module fast_loc
!!===============================================================================84
!!
!!  fast localization series fuctions.
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2018.03.13               hanrui                   original code
!!        2018.04.01               Hanrui             NV2.0 Final version
!!                                                  
!!===============================================================================84

    use common,             only : nx,ny,nbox
    use common,             only : rx,ry
    use inputoption,        only : Loc_uv
    use common,             only : MLonC,MLatC
    use interpolation,      only : SpheDistance
    use Mathmatic,          only : SVD_Lapack_r8
    use readwrite_nc,       only : Write_CO2_nc
    !   reshape
    use reshapeTT,          only : reshape2Dto1D
    use reshapeTT,          only : reshape1Dto2D
    use reshapeTT,          only : reshape3Dto1D
    use reshapeTT,          only : reshape1Dto3D
    !   ocogst_operator
    use operator_ocogst,    only : locate_startpoint
    use operator_ocogst,    only : locate_Boxaxis
    use operator_ocogst,    only : locate_weight

    implicit None
    !--------------------------------------------------------------------------
    PUBLIC                  ::  RouO_Operator
    PUBLIC                  ::  FormRxy
    !   gen Cx,RouX
    PUBLIC                  ::  FormCx
    PUBLIC                  ::  FormRx
    !   gen Cy,RouY
    PUBLIC                  ::  FormCy
    PUBLIC                  ::  FormRy
    !   caculate rou
    PUBLIC                  ::  Rou

    !--------------------------------------------------------------------------
    contains

    subroutine RouO_Operator( nobs, lon, lat, RouXY, RouObs)
    !!=======================================================================76
    !!
    !!      from RouXYZ to RouObs
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            integer,intent(in)          ::  nobs
            real(8),intent(in)          ::  lon(nobs)
            real(8),intent(in)          ::  lat(nobs)
            real(8),intent(in)          ::  RouXY(nx,ny,rx*ry)
            real(8),intent(out)         ::  RouObs(nobs,rx*ry)

            integer                     ::  ii, jj, kk, ll

            integer                     ::  stpt_lon, stpt_lat
            real(8)                     ::  Roubox(nbox,nbox,rx*ry)
            real(8)                     ::  Box_lon(nbox), Point_lon
            real(8)                     ::  Box_lat(nbox), Point_lat
            real(8)                     ::  Box_wght( nbox, nbox)
            !----------------------------------------------------------
            RouObs = 0.0d0

            do ll=1,nobs
            !   locate start point (checked)
            call locate_startpoint( real(lon(ll)), real(lat(ll)), stpt_lon, stpt_lat )
            call locate_Boxaxis( real(lon(ll)), real(lat(ll)), &
                &Box_lon, Box_lat, Point_lon, Point_lat )

            if( stpt_lon.lt.nx .and. stpt_lon.ge.1 ) then
                Roubox(:,:,:)  = RouXY( stpt_lon:stpt_lon+1, stpt_lat:stpt_lat+1, : )
            elseif( stpt_lon .eq. nx ) then
                Roubox(1,:,:)  = RouXY( nx, stpt_lat:stpt_lat+1, : )
                Roubox(2,:,:)  = RouXY( 1,  stpt_lat:stpt_lat+1, : )
            endif

            call locate_weight(Box_lon, Box_lat, Point_lon, Point_lat, Box_wght)

            do ii=1,rx*ry
            RouObs(ll,ii) = sum( Box_wght(:,:) * Roubox(:,:,ii) )
            enddo

            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine RouO_Operator

    subroutine FormRxy( RouX, RouY, RouXY )
    !!=======================================================================76
    !!
    !!      to generate Rxy ( nx*ny, rx*ry ).  org.ZHQ. FormShensiRxy
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real(8),intent(in)      ::  RouX(nx,rx,ny)
            real(8),intent(in)      ::  RouY(ny,ry)
            real(8),intent(out)     ::  RouXY(nx*ny,rx*ry)


            real(8)                 ::  SPx(nx*ny,rx)
            real(8)                 ::  PRy(nx*ny,ry)

            integer                 ::  ii, jj, kk

            !----------------------------------------------------------

            do ii=1,rx
            call reshape2Dto1D(nx,ny,RouX(:,ii,:),SPx(:,ii))
            enddo

            do jj = 1,ny
            do ii = 1,nx
                PRy((jj-1)*nx+ii,:) = RouY(jj,:)
            enddo
            enddo

            do jj = 1,ry
            do ii = 1,rx
              RouXY( :, (jj-1)*rx + ii ) = PRy( :, jj) * SPx( :, ii)
            enddo
            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine FormRxy


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine FormRx( RouX )
    !!=======================================================================76
    !!
    !!      to generate Rx ( ny, nx, rx ).  org.ZHQ. mod_x
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real(8),intent(out)     ::  RouX(nx,rx,ny)
            !
            real(8)                 ::  lat
            real(8)                 ::  Cx(nx,nx)
            real(8)                 ::  Ux(nx,nx)
            real(8)                 ::  UxT(nx,nx)
            real(8)                 ::  SigMa(nx,nx)
            !
            real(8)                 ::  tmp_lam01(nx,rx)
            real(8)                 ::  tmp_lam02(rx,rx)
            real(8)                 ::  lam(rx,rx)
            !
            integer                 ::  ii, jj, kk, ll

            !----------------------------------------------------------
            call FormCx( MLatC( floor(( ny ) / 2.0 ) + 1), Cx )

            !   use Ux(:,1:rx), UxT(1:rx,:)
            call SVD_Lapack_r8( nx, nx, Cx, Ux, SigMa, UxT )

            do ii = 1,ny

                lat = MLatC(ii)

                call FormCx( lat, Cx )

                tmp_lam01(:,:) = matmul( Cx, Ux(:,1:rx) )

                tmp_lam02(:,:) = matmul( UxT(1:rx,:), tmp_lam01 )

                lam = 0.0d0

                do jj=1,rx
                lam(jj,jj) = sqrt( tmp_lam02(jj,jj) )
                enddo

                RouX( :, :, ii) = matmul( Ux(:,1:rx), lam ) 

            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine FormRx

    subroutine FormCx( lat, Cx )
    !!=======================================================================76
    !!
    !!      to generate Cx ( nx, nx ) of lat.  org.ZHQ. FormCx
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in and out
            real(8),intent(in)  ::  lat
            real(8),intent(out) ::  Cx( nx, nx )

            real(8)             ::  dis

            integer             ::  ii, jj, kk

            !----------------------------------------------------------

            do ii = 1,nx
            do jj = 1,nx

                dis = SpheDistance( MLonC(ii), lat, MLonC(jj), lat)
                dis = dis/Loc_uv
                Cx(ii,jj) = Rou(dis)

            enddo
            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine FormCx

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine FormRy( RouY )
    !!=======================================================================76
    !!
    !!      to generate RouY ( ny, ry ).  org.ZHQ. mod_y
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in out
            real(8),intent(out) ::  RouY(ny,ry)
            !
            real(8)             ::  Cy(ny,ny)
            real(8)             ::  Uy(ny,ny)
            real(8)             ::  UyT(ny,ny)
            real(8)             ::  SigMa(ny,ny)

            integer             ::  ii, jj, kk
            !----------------------------------------------------------

            call FormCy( Cy )

            call SVD_Lapack_r8( ny, ny, Cy, Uy, SigMa, UyT)

            do ii = 1 , ry
            RouY(:,ii) = Uy(:,ii) * sqrt( SigMa( ii, ii ) )
            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine FormRy


    subroutine FormCy( Cy )
    !!=======================================================================76
    !!
    !!      to generate Cy ( ny, ny ).  org.ZHQ. FormCy
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   in and out
            real(8),intent(out) ::  Cy(ny,ny)

            real(8)             ::  dis

            integer             ::  ii, jj, kk
            !----------------------------------------------------------

            do ii = 1,ny
            do jj = 1,ny

                dis = SpheDistance( 0.0d0, MLatC(ii), 0.0d0, MLatC(jj))
                dis = dis/Loc_uv
                Cy(ii,jj) = Rou(dis)

            enddo
            enddo

            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine FormCy

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Rou(dis)
    !!=======================================================================76
    !!
    !!      Rou(dis)
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input and output 
            real(8),intent(in)      ::  dis
            real(8)                 ::  Rou

            !----------------------------------------------------------
            !   parameter
            integer                 ::  i,      j,      k
            integer                 ::  ii,     jj,     kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            if( (dis.ge.0.0d0) .and. (dis.lt.1.0d0) ) then 

                Rou = 1.0d0-(dis**5/4.0d0)+(dis**4/2.0d0)
                Rou = Rou+(5*dis**3/8.0d0)-(5*dis**2/3.0d0)

            else if( (dis.ge.1.0d0) .and. (dis.lt.2.0d0) ) then

                Rou = 4.0d0+(dis**5/12.0d0)-(dis**4/2.0d0)
                Rou = Rou+(5*dis**3/8.0d0)+(5*dis**2/3.0d0)
                Rou = Rou-5*dis-2.0d0/(dis*3.0d0)

            else
            
                Rou = 0.0d0
            
            endif

            return

            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end function Rou


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
End Module fast_loc
