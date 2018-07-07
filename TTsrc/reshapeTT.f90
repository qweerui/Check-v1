MODULE ReshapeTT
!!=================================================================================
!!
!!    Input parameters
!!
!!    Record of revision:
!!           Date:              Programmer:              Description of change:
!!     ------------           --------------            ------------------------
!!     2017.10.29                   Hanrui                 NV2.0 original code
!!     2018.03.14                   Hanrui                   NV2.0 MDF for loc
!!     2018.03.28                   Hanrui                 NV2.0 Final version
!!
!!=================================================================================
    IMPLICIT NONE

    !--------------------------------------------------------------------------
    PUBLIC          ::  reshape2Dto1D
    PUBLIC          ::  reshape1Dto2D
    PUBLIC          ::  reshape3Dto1D
    PUBLIC          ::  reshape1Dto3D

    !--------------------------------------------------------------------------
    Contains
    
    subroutine reshape2Dto1D(xx,yy,Matrix,Matrix1D)
    !========================================================================80
    !
    !   Reshape 2dimensional matrix Matrix(xx,yy) to 1dimensional Matrix1D(xx*yy)
    !
    !                                                   Rui 2016/07/12
    !========================================================================80
        implicit none
        !------------------------------------------------------------------
            integer,intent(in)      ::  xx
            integer,intent(in)      ::  yy
            real(8),intent(in)      ::  Matrix(xx,yy)
            real(8),intent(out)     ::  Matrix1D(xx*yy)
            integer                 ::  i,  j

        !------------------------------------------------------------------
            ! do i=1,xx
            ! do j=1,yy
            !     Matrix1D( yy*(i-1)+j ) = Matrix(i,j)
            ! enddo
            ! enddo

            !for loc
            do j=1,yy
                Matrix1D( xx*(j-1)+1 : xx*j ) = Matrix( : , j )
            enddo

        !------------------------------------------------------------------

    !========================================================================80
    end subroutine reshape2Dto1D

    subroutine reshape1Dto2D(xx,yy,Matrix1D,Matrix)
    !========================================================================80
    !
    !   Reshape 1dimensional Matrix1D(xx*yy) to 2dimensional matrix Matrix(xx,yy)
    !
    !                                                   Rui 2016/07/12
    !========================================================================80
        implicit none
        !------------------------------------------------------------------
            integer,intent(in)      ::  xx
            integer,intent(in)      ::  yy
            real(8),intent(in)      ::  Matrix1D(xx*yy)
            real(8),intent(out)     ::  Matrix(xx,yy)
            integer                 ::  i,  j

        !------------------------------------------------------------------
            ! do i=1,xx
            ! do j=1,yy
            !     Matrix(i,j) = Matrix1D( yy*(i-1)+j )
            ! enddo
            ! enddo

            !for loc
            do j=1,yy
                Matrix( : , j ) = Matrix1D( xx*(j-1)+1 : xx*j )
            enddo

        !------------------------------------------------------------------

    !========================================================================80
    end subroutine reshape1Dto2D

    subroutine reshape3Dto1D(xx,yy,zz,Matrix,Matrix1D)
    !========================================================================80
    !
    !   Reshape 3dimensional matrix Matrix(xx,yy,zz) to 
    !   1dimensional Matrix1D(xx*yy*zz)
    !
    !                                                   Rui 2016/07/12
    !========================================================================80
        implicit none
        !------------------------------------------------------------------
            integer,intent(in)      ::  xx
            integer,intent(in)      ::  yy
            integer,intent(in)      ::  zz
            real(8),intent(in)      ::  Matrix(xx,yy,zz)
            real(8),intent(out)     ::  Matrix1D(xx*yy*zz)
            real(8)                 ::  Matrix2D(xx*yy,zz)
            integer                 ::  i,  j,  k

        !------------------------------------------------------------------
            ! do i=1,xx
            ! do j=1,yy
            ! do k=1,zz
            !     Matrix1D( zz*((yy*(i-1)+j)-1)+k ) = Matrix(i,j,k)
            ! enddo
            ! enddo
            ! enddo

            ! for loc
            do k=1,zz
            call reshape2Dto1D( xx, yy, Matrix(:,:,k), Matrix2D(:,k))
            enddo
            call reshape2Dto1D( xx*yy, zz, Matrix2D, Matrix1D)

        !------------------------------------------------------------------

    !========================================================================80
    end subroutine reshape3Dto1D

    subroutine reshape1Dto3D(xx,yy,zz,Matrix1D,Matrix)
    !========================================================================80
    !
    !   Reshape 1dimensional Matrix1D(xx*yy*zz) to 
    !   3dimensional matrix Matrix(xx,yy,zz)
    !
    !                                                   Rui 2016/07/12
    !========================================================================80
        implicit none
        !------------------------------------------------------------------
            integer,intent(in)      ::  xx
            integer,intent(in)      ::  yy
            integer,intent(in)      ::  zz
            real(8),intent(in)      ::  Matrix1D(xx*yy*zz)
            real(8),intent(out)     ::  Matrix(xx,yy,zz)
            real(8)                 ::  Matrix2D(xx*yy,zz)
            integer                 ::  i,  j,  k

        !------------------------------------------------------------------
            ! do i=1,xx
            ! do j=1,yy
            ! do k=1,zz
            !     Matrix(i,j,k) = Matrix1D( zz*((yy*(i-1)+j)-1)+k )
            ! enddo
            ! enddo
            ! enddo

            ! for loc
            call reshape1Dto2D( xx*yy, zz, Matrix1D, Matrix2D)
            do k=1,zz
            call reshape2Dto1D( xx, yy, Matrix2D(:,k), Matrix(:,:,k))
            enddo

        !------------------------------------------------------------------

    !========================================================================80
    end subroutine reshape1Dto3D

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!==================================================================================
end MODULE ReshapeTT
