Module PODCORE
!!===============================================================================84
!!  contains:
!!
!!      subroutine PODEn4DVarCORE
!!
!!      subroutine EnsembleUpdate
!!
!!      subroutine localization_H
!!
!!      subroutine localization_V
!!
!!      function Rou(dis)
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2016.09.04               hanrui                   original code
!!
!!        2016.12.24               hanrui              MDF. PODEn4DVarCORE
!!                    change variable AA(ndim,nrens) and obsAA(nobs,nrens)
!!                    from intent(inout) to intent(in) add new variable 
!!                    AAptb(ndim,nrens) and ObsAAptb(nobs,nrens)
!!
!!        2016.12.24               hanrui              ORG. EnsembleUpdate
!!        2016.07.21               hanrui                  ORG. PODCOREpre
!!                                                  
!!===============================================================================84

    implicit None
    !--------------------------------------------------------------------------
    PUBLIC          ::  PODEn4DVarCORE
                    !   PODEn4DVar Core
                    !   use Mathmatic,          only : SVD_Lapack_r8,inv_r8
                    !   use reshapeTT,          only : reshape3Dto1D
    PUBLIC          ::  EnsembleUpdate
                    !   Ensemble Update
                    !   use Mathmatic,          only : SVD_Lapack_r8,inv_r8
                    !   use reshapeTT,          only : reshape3Dto1D
    PUBLIC          ::  localization_H
                    !   horizontal localization
                    !   use Common,             only : nx,ny,MLonC,MLatC
                    !   use interpolation,      only : SpheDistance
                    !   use InputOption,        only : Loc_uv
    PUBLIC          ::  localization_V
                    !   vertical localization
                    !   use Common,             only : nz
                    !   use Common,             only : LevAlt
                    !   use InputOption,        only : Loc_uv

    !--------------------------------------------------------------------------
    contains

    subroutine PODEn4DVarCORE(     nx , ny , nz ,  nrens ,   nobs  ,&
                            &     AA  ,    X_b  ,  ObsAA ,  ObsXb  ,&
                            & ObsVal  , ObsErr  , ObsAlt , ObsLat  ,&
                            & ObsLon  ,    X_a  , X_aErr )
    !!=======================================================================76
    !!
    !!  PODEn4DVar Core
    !!
    !!=======================================================================76
        use Mathmatic,          only : SVD_Lapack_r8
        use Mathmatic,          only : inv_r8
        use reshapeTT,          only : reshape3Dto1D
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input
            integer,intent(in)      ::  nx
            integer,intent(in)      ::  ny
            integer,intent(in)      ::  nz
            integer,intent(in)      ::  nrens
            integer,intent(in)      ::  nobs
            real(8),intent(in)      ::  AA(nx*ny*nz,nrens)
            real(8),intent(in)      ::  X_b(nx*ny*nz)
            real(8),intent(in)      ::  ObsAA(nobs,nrens)
            real(8),intent(in)      ::  ObsXb(nobs)
            real(8),intent(in)      ::  ObsVal(nobs)
            real(8),intent(in)      ::  ObsErr(nobs)
            real(8),intent(in)      ::  ObsAlt(nobs)
            real(8),intent(in)      ::  ObsLat(nobs)
            real(8),intent(in)      ::  Obslon(nobs)

            !   ptb
            real(8)                 ::  AAptb(nx*ny*nz,nrens)
            real(8)                 ::  ObsAAptb(nobs,nrens)

            
            !   inner
            real(8)                 ::  ObsXbErr(nobs)

            !   output
            real(8),intent(out)     ::  X_a(nx*ny*nz)
            real(8),intent(out)     ::  X_aErr(nx*ny*nz)

            !   parameter
            integer                 ::  ii, jj, kk
            integer                 ::  i,  j,  k

            !   singliar value decomposition
            real(8)                 ::  TT(nrens,nrens)
            real(8)                 ::  VV(nrens,nrens)
            real(8)                 ::  SigMa(nrens,nrens)
            real(8)                 ::  VVT(nrens,nrens)

            !   PODEn4DVar parameter
            integer                 ::  Tnrens
            real(8),allocatable     ::  Phiy(:,:)       !nobs,Tnrens
            real(8),allocatable     ::  PhiyT(:,:)      !Tnrens,nobs
            real(8),allocatable     ::  Psi_A(:,:)      !Tnrens,Tnrens
            real(8),allocatable     ::  PhiywaveT(:,:)  !Tnrens,nobs
            real(8),allocatable     ::  Phix(:,:)       !nx*ny*nz,Tnrens
            real(8),allocatable     ::  PhixT(:,:)      !Tnrens,nx*ny*nz
            real(8)                 ::  PhixPhiywaveT(nx*ny*nz)

            !   localization
            real(8)                 ::  Rou_z(nz)
            real(8)                 ::  Rou_xy(nx,ny)
            real(8)                 ::  Rou_3D(nx,ny,nz)
            real(8)                 ::  Rou_xyz(nx*ny*nz)

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   samples perturbation
            !           y' = Samples_y - Xb_y
            !           y'obs = Obs - ObsXb
            do ii=1,nobs
            do jj=1,nrens
                ObsAAptb(ii,jj) = ObsAA(ii,jj) - ObsXb(ii)
            enddo
                ObsXbErr(ii) = ObsVal(ii) - ObsXb(ii)
            enddo

            write(*,*) " Samples perturbation finished! "
            write(*,*) ""

            !----------------------------------------------------------
            !   Single Value Decomposition
            !       TT = y'T * y' = VV SigMa VV T
            TT = matmul(transpose(ObsAAptb),ObsAAptb)
            call SVD_Lapack_r8(nrens,nrens,TT,VV,SigMa,VVT)

            write(*,*) " y'T * y' = VV * SigMa * VVT finished! "
            write(*,*) ""

            !   allocate
            Tnrens = nrens
            allocate(Phiy(nobs,Tnrens),PhiyT(Tnrens,nobs))
            allocate(Psi_A(Tnrens,Tnrens),PhiywaveT(Tnrens,nobs))
            allocate(Phix(nx*ny*nz,Tnrens),PhixT(Tnrens,nx*ny*nz))

            !----------------------------------------------------------
            !       Phiy = y' * VV
            Phiy = matmul(ObsAAptb,VV(:,1:Tnrens))
            PhiyT = transpose(Phiy)

            write(*,*) " Phiy = y' * VV finished"
            write(*,*) ""

            !----------------------------------------------------------
            !       Psi_A = [(Tnrens-1)I + PhiyT * R-1 * Rhiy] -1
            !       PhiywaveT = Psi_A * PhiyT * R-1
            Psi_A = 0.0d0
            do jj=1,Tnrens
                Psi_A(jj,jj) = (Tnrens-1)*1.0d0
            do ii=1,nobs
                PhiyT(jj,ii) = PhiyT(jj,ii) / (ObsErr(ii)*ObsErr(ii))
            enddo
            enddo
            Psi_A = Psi_A + matmul(PhiyT,Phiy)
            Psi_A = inv_r8(Psi_A)
            PhiywaveT = matmul(Psi_A,PhiyT)

            write(*,*) " Psi_A = [(Tnrens-1)I + PhiyT * R-1 * Rhiy] -1 finished"
            write(*,*) ""
            write(*,*) " PhiywaveT = Psi_A * PhiyT * R-1 finished"
            write(*,*) ""

            !----------------------------------------------------------
            !       AAptb = Simple - background
            !       Phix = AAptb * VV
            !       X_aErr(i) = Phix * Psi_A * PhixT (i,i)
            !       X_a
            do jj = 1,nrens
            do kk = 1,nx*ny*nz
            AAptb(kk,jj) = AA(kk,jj) - X_b(kk)
            enddo
            enddo

            Phix  = matmul(AAptb,VV(:,1:Tnrens))
            PhixT = transpose(Phix)
            PhixT = matmul(Psi_A,PhixT)

            do kk=1,nx*ny*nz
                X_aErr(kk) = sqrt(sum(Phix(kk,:)*PhixT(:,kk)))
            enddo
            X_a = 0.0d0

            do ii=1,nobs
                PhixPhiywaveT = 0.0d0
                write(*,*) "Obs cycle:",ii
                ! - - - - - - - - - - - - - - - - - - - - - - - - -
                !   horizontal localization
                call localization_H(ObsLon(ii),ObsLat(ii),Rou_xy)
                !   vertical localization
                !call localization_V(nz,ObsAlt(ii),Rou_z)
                do jj=1,nz
                Rou_3D(:,:,jj) = Rou_xy(:,:)
                !Rou_3D(:,:,jj) = Rou_3D(:,:,jj) * Rou_z(jj)
                enddo
                call reshape3Dto1D(nx,ny,nz,Rou_3D,Rou_xyz)
                ! - - - - - - - - - - - - - - - - - - - - - - - - -
                do kk=1,nx*ny*nz
                PhixPhiywaveT(kk) = sum(phix(kk,:)*PhiywaveT(:,ii))
                X_a(kk) = X_a(kk) + Rou_xyz(kk)*PhixPhiywaveT(kk)*ObsXbErr(ii)
                enddo
            enddo
            X_a = X_a + X_b

            !----------------------------------------------------------
            !   clear
            deallocate(Phiy,PhiyT)
            deallocate(Psi_A,PhiywaveT)
            deallocate(Phix,PhixT)
            return

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine PODEn4DVarCORE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine EnsembleUpdate(      nx ,  ny ,  nz ,  nrens  ,  nobs  ,   AA ,&
                            & Inflation, X_b , ObsAA , ObsXb ,ObsErr  ,   X_a)
    !!=======================================================================76
    !!
    !!  EnsembleUpdate
    !!      samename:    AA(input)      ->      AA(output)
    !!
    !!=======================================================================76
        use Mathmatic,          only : SVD_Lapack_r8
        use Mathmatic,          only : inv_r8
        use reshapeTT,          only : reshape3Dto1D
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input
            integer,intent(in)      ::  nx
            integer,intent(in)      ::  ny
            integer,intent(in)      ::  nz
            integer,intent(in)      ::  nrens
            integer,intent(in)      ::  nobs
            real(8),intent(inout)   ::  AA(nx*ny*nz,nrens)
            real(8),intent(in)      ::  Inflation
            real(8),intent(in)      ::  X_b(nx*ny*nz)
            real(8),intent(in)      ::  ObsAA(nobs,nrens)
            real(8),intent(in)      ::  ObsXb(nobs)
            real(8),intent(in)      ::  ObsErr(nobs)
            real(8),intent(in)      ::  X_a(nx*ny*nz)

            !   ptb
            real(8)                 ::  AAptb(nx*ny*nz,nrens)
            real(8)                 ::  ObsAAptb(nobs,nrens)
            
            !   parameter
            integer                 ::  ii, jj, kk
            integer                 ::  i,  j,  k

            !   singliar value decomposition
            real(8)                 ::  TT(nrens,nrens)
            real(8)                 ::  VV(nrens,nrens)
            real(8)                 ::  SigMa(nrens,nrens)
            real(8)                 ::  VVT(nrens,nrens)

            !   Psi_A parameter
            real(8)                 ::  Phiy(nobs,nrens)
            real(8)                 ::  PhiyT(nrens,nobs)
            real(8)                 ::  Psi_A(nrens,nrens)

            !   output
            real(8)                 ::  AAnew_ptb(nx*ny*nz,nrens)

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   samples perturbation
            !           y' = Samples_y - Xb_y
            !           y'obs = Obs - ObsXb
            do ii=1,nobs
            do jj=1,nrens
                ObsAAptb(ii,jj) = ObsAA(ii,jj) - ObsXb(ii)
            enddo
            enddo

            write(*,*) " Samples perturbation finished! "
            write(*,*) ""

            !----------------------------------------------------------
            !   Single Value Decomposition
            !       TT = y'T * y' = VV SigMa VV T
            TT = matmul(transpose(ObsAAptb),ObsAAptb)
            call SVD_Lapack_r8(nrens,nrens,TT,VV,SigMa,VVT)

            write(*,*) " y'T * y' = VV * SigMa * VVT finished! "
            write(*,*) ""

            !----------------------------------------------------------
            !       Phiy = y' * VV
            Phiy  = matmul(ObsAAptb,VV(:,1:nrens))
            PhiyT = transpose(Phiy)

            write(*,*) " Phiy = y' * VV finished"
            write(*,*) ""

            !----------------------------------------------------------
            !       Psi_A = [(nrens-1)I + PhiyT * R-1 * Rhiy] -1
            Psi_A = 0.0d0
            do jj=1,nrens
                Psi_A(jj,jj) = (nrens-1)*1.0d0
            do ii=1,nobs
                PhiyT(jj,ii) = PhiyT(jj,ii) / (ObsErr(ii)*ObsErr(ii))
            enddo
            enddo
            Psi_A = Psi_A + matmul(PhiyT,Phiy)
            Psi_A = inv_r8(Psi_A)

            write(*,*) " Psi_A = [(nrens-1)I + PhiyT * R-1 * Rhiy] -1 finished"
            write(*,*) ""

            !----------------------------------------------------------
            !   AAptb = Simple - background
            do jj = 1,nrens
            do kk = 1,nx*ny*nz
            AAptb(kk,jj) = AA(kk,jj) - X_b(kk)
            enddo
            enddo

            !----------------------------------------------------------
            !   AAnew_ptb = AAptb * ( (nrens-1)Psi_A ) ^ 1/2
            !   AA = Inflation * AAnew_ptb + X_a
            Psi_A = Psi_A * (nrens-1)
            call SVD_Lapack_r8(nrens,nrens,Psi_A,VV,SigMa,VVT)
            do i=1,nrens
            SigMa(i,i) = sqrt(SigMa(i,i))
            enddo
            Psi_A = matmul(VV,SigMa)
            AAnew_ptb = matmul(AAptb,Psi_A)
            do i=1,nrens
            AA(:,i) = Inflation*AAnew_ptb(:,i)+X_a(:)
            enddo

            write(*,*) " AAnew_ptb = AAptb * ( (nrens-1)Psi_A ) ^ 1/2"
            write(*,*) " AA = Inflation * AAnew_ptb + X_a"
            write(*,*) ""

            !----------------------------------------------------------
            return

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine EnsembleUpdate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine localization_H(longitude,latitude,Rou_xy)
    !!=======================================================================76
    !!
    !!  caculate loclization matrix of obs(jj)
    !!
    !!=======================================================================76
        use Common,             only : nx,ny
        use Common,             only : MLonC,MLatC
        use interpolation,      only : SpheDistance
        use InputOption,        only : Loc_uv
        use fast_loc,           only : Rou
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input and output 
            real(8),intent(in)      ::  longitude
            real(8),intent(in)      ::  latitude
            real(8),intent(out)     ::  Rou_xy(nx,ny)

            !----------------------------------------------------------
            !   parameter
            real(8)                 ::  dis
            integer                 ::  ii,     jj,     kk

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            Rou_xy = 0.0d0
            do ii=1,nx
            do jj=1,ny
            dis = SpheDistance(longitude,latitude,MLonC(ii),MLatC(jj))
            dis = dis/Loc_uv
            Rou_xy(ii,jj) = Rou(dis)
            enddo
            enddo
            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine localization_H

    subroutine localization_V(nn,Altitude,Rou_z)
    !!=======================================================================76
    !!
    !!  wrong!
    !!      not allowed to use!
    !!
    !!  caculate loclization matrix of obs(jj)
    !!
    !!=======================================================================76
        use Common,             only : nz
        use Common,             only : LevAlt
        use InputOption,        only : Loc_uv
        use fast_loc,           only : Rou
        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input and output 
            integer,intent(in)      ::  nn
            real(8),intent(in)      ::  Altitude
            real(8),intent(out)     ::  Rou_z(nz)

            !----------------------------------------------------------
            !   parameter
            real(8)                 ::  dis
            integer                 ::  ii

            !----------------------------------------------------------

        !------------------------------------------------------------------

            stop "wrong localization_V"

            !----------------------------------------------------------
            !   for flx only
            if(nn.ne.nz) then
                Rou_z(:) = 1.0d0
                return
            endif

            !   for gosat obs without altitude
            if(Altitude.eq.-999.0d0) then
                Rou_z(:) = 1.0d0
                return
            endif

            Rou_z = 0.0d0
            do ii=1,nz
            dis = abs(altitude - LevAlt(ii))/Loc_uv
            Rou_z(ii) = Rou(dis)
            enddo
            return
            !----------------------------------------------------------

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine localization_V

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
End Module PODCORE
