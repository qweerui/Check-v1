Module NLSCORE
!!===============================================================================84
!!  contains:
!!
!!      subroutine NLS5_4DVarCORE
!!
!!      Record of revision:
!!              Date:          Programmer:          Description of change:
!!        ------------       --------------        ------------------------
!!        2018.03.12              hanrui                   original code
!!        2018.04.01              Hanrui             NV2.0 Final version
!!                                                  
!!===============================================================================84

    use Mathmatic,          only : SVD_Lapack_r8
    use Mathmatic,          only : inv_r8
    use reshapeTT,          only : reshape3Dto1D
    use InputOption,        only : DA_obs_bmc
    use common,             only : rx,ry,rz
    use fast_loc

    implicit None
    !--------------------------------------------------------------------------
    !   NLS5 4DVar Core
    PUBLIC                  ::  NLS5_4DVarCORE
    !   NLS5 Ensemble Update
    PUBLIC                  ::  NLS5_EsbUpdate

    contains

    subroutine NLS5_4DVarCORE( nx, ny, nz, nrens, nobs, DA_obs, AA, X_b, X_a, X_aErr)
    !!=======================================================================76
    !!
    !!  NLS5 4DVar Core
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input
            integer,intent(in)          ::  nx
            integer,intent(in)          ::  ny
            integer,intent(in)          ::  nz
            integer,intent(in)          ::  nrens
            integer,intent(in)          ::  nobs
            type(DA_obs_bmc),intent(in) ::  DA_obs(nobs)
            real(8),intent(in)          ::  AA(nx*ny*nz,nrens)
            real(8),intent(in)          ::  X_b(nx*ny*nz)

            !   output
            real(8),intent(out)         ::  X_a(nx*ny*nz)
            real(8),intent(out)         ::  X_aErr(nx*ny*nz)

            !   ptb
            real(8)                     ::  AAptb(nx*ny*nz,nrens)
            real(8)                     ::  ObsAAptb(nobs,nrens)
            
            !   inner
            real(8)                     ::  ObsXbErr(nobs)

            !   parameter
            integer                     ::  ii, jj, kk, ll, mm
            integer                     ::  i,  j,  k,  l,  m

            !   singliar value decomposition
            real(8)                     ::  TT(nrens,nrens)
            real(8)                     ::  VV(nrens,nrens)
            real(8)                     ::  SigMa(nrens,nrens)
            real(8)                     ::  VVT(nrens,nrens)

            !   PODEn4DVar parameter
            integer                     ::  Tnrens
            real(8),allocatable         ::  Phiy(:,:)       !nobs,Tnrens
            real(8),allocatable         ::  PhiyT(:,:)      !Tnrens,nobs
            real(8),allocatable         ::  Psi_A(:,:)      !Tnrens,Tnrens
            real(8),allocatable         ::  Phiywave(:,:)   !nobs,Tnrens
            real(8),allocatable         ::  PhiywaveT(:,:)  !Tnrens,nobs
            real(8),allocatable         ::  Phix(:,:)       !nx*ny*nz,Tnrens
            real(8),allocatable         ::  PhixT(:,:)      !Tnrens,nx*ny*nz
            real(8),allocatable         ::  PhiywaveobsTmp(:)

            !   localization

            real(8),allocatable         ::  PhixRou_Row(:)
            real(8)                     ::  PhiywaveRou_Col(nobs)
            real(8)                     ::  RouXYZ_Row(rx*ry*rz)
            real(8),allocatable         ::  PhiywaveYptb_Row(:)

            !   localization
            real(8)                     ::  RouX(nx,rx,ny)
            real(8)                     ::  RouY(ny,ry)
            real(8)                     ::  RouXY(nx*ny,rx*ry)
            real(8)                     ::  RouXY_3Dtmp(nx,ny,rx*ry)
            !
            real(8)                     ::  RouZ(nz,rz)
            real(8)                     ::  RouObs(nobs,rx*ry*rz)

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   samples perturbation
            !           y' = Samples_y - Xb_y
            !           y'obs = Obs - ObsXb
            do ii=1,nobs
            do jj=1,nrens
                ObsAAptb(ii,jj) = DA_obs(ii)%obs_AA(jj) - DA_obs(ii)%obs_Xb
            enddo
                ObsXbErr(ii) = DA_obs(ii)%obs - DA_obs(ii)%obs_Xb
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
            allocate(Psi_A(Tnrens,Tnrens))
            allocate(Phiywave(nobs,Tnrens),PhiywaveT(Tnrens,nobs))
            allocate(Phix(nx*ny*nz,Tnrens),PhixT(Tnrens,nx*ny*nz))

            allocate(PhixRou_Row(Tnrens*rx*ry*rz))
            allocate(PhiywaveYptb_Row(rx*ry*rz*Tnrens))

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
                PhiyT(jj,ii) = PhiyT(jj,ii) / ( DA_obs(ii)%obs_err * DA_obs(ii)%obs_err )
            enddo
            enddo
            Psi_A = Psi_A + matmul(PhiyT,Phiy)
            Psi_A = inv_r8(Psi_A)
            PhiywaveT = matmul(Psi_A,PhiyT)
            Phiywave  = transpose(PhiywaveT)

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

            ! !----------------------------------------------------------
            ! !   no localization
            ! ! do ii=1,nobs
            ! ! do kk=1,nx*ny*nz
            ! ! PhiywaveobsTmp = sum(phix(kk,:)*PhiywaveT(:,ii))
            ! ! X_a(kk) = X_a(kk) + PhiywaveobsTmp*ObsXbErr(ii)
            ! ! enddo
            ! ! enddo

            ! allocate(PhiywaveobsTmp(Tnrens))
            ! do ii=1,Tnrens
            ! PhiywaveobsTmp(ii) = sum(PhiywaveT(ii,:)*ObsXbErr(:))
            ! enddo
            ! do kk=1,nx*ny*nz
            ! X_a(kk) =  sum(Phix(kk,:)*PhiywaveobsTmp(:))
            ! enddo
            ! deallocate(PhiywaveobsTmp)
            ! X_a = X_a + X_b

            ! deallocate(Phiy,PhiyT)
            ! deallocate(Psi_A,PhiywaveT)
            ! deallocate(Phix,PhixT)
            ! return

            !----------------------------------------------------------
            !   localization
            write(*,*) "Generate RouX:"
            call FormRx( RouX )
            write(*,*) "Generate RouY:"
            call FormRy( RouY )
            write(*,*) "Generate RouXY:"
            call FormRxy ( RouX, RouY, RouXY )
            write(*,*) "Generate RouZ:"
            RouZ(:,:) = 1.0d0/sqrt(rz*1.0d0)

            write(*,*) "Generate RouObs:"
            do jj=1,rz
            write(*,*) "Finish ",   jj, "th"
            do ii=1,rx*ry
            call reshape1Dto2D( nx, ny, RouXY(:,ii), RouXY_3Dtmp(:,:,ii) )
            enddo
            call RouO_Operator( nobs, DA_obs(:)%lon, DA_obs(:)%lat, &
                &RouXY_3Dtmp(:,:,:), RouObs( :, (jj-1)*rx*ry+1 : jj*rx*ry ) )
            enddo
            write(*,*) "Generate finish."

            !   fast localization schame:
            !       ( X ) : Matrix that cannot save.
            !
            !   RouXY,          RouZ,       RouObs,             RouXYZ
            !   nx*ny x rx*ry,  nz x rz,    nobs x rx*ry*rz,    nx*ny*nz x rx*ry*rz ( X )
            !
            !                 nobs x rx*ry*rz
            !   PhiywaveRou = RouObs    <e>    Phiywave
            !   nobs x Tnrens*rx*ry*rz ( X )   nobs x Tnrens
            !
            !             nx*ny*nz x rx*ry*rz
            !   PhixRou = RouXYZ       <e>       Phix
            !   nx*ny*nz x Tnrens*rx*ry*rz ( X ) nx*ny*nz x Tnrens
            !
            !   PhixPhiywave = mutmal ( PhixRou , transpose( PhiywaveRou ) )
            !   nx*ny*nz x nobs ( X )
            !
            !   X_a = mutmal ( PhixPhiywave , ObsXbErr )
            !   nx*ny*nz

            do mm = 1,Tnrens
            do ll = 1,rx*ry*rz

                PhiywaveRou_Col(:) = RouObs(:,ll) * Phiywave(:,mm)

                jj = ( mm - 1 ) * rx*ry*rz + ll
                
                PhiywaveYptb_Row(jj) = sum( PhiywaveRou_Col(:) * ObsXbErr(:) )

            enddo
            enddo

            do mm = 1,nz
            do ll = 1,nx*ny

                do kk=1,rz
                RouXYZ_Row( (kk-1)*rx*ry+1 : kk*rx*ry ) = RouXY( ll, : ) * RouZ( mm, kk )
                enddo

                ii = ( mm - 1 ) * nx*ny + ll

                do kk=1,Tnrens
                PhixRou_Row( (kk-1)*rx*ry*rz+1 : kk*rx*ry*rz ) = RouXYZ_Row(:) * Phix(ii,kk)
                enddo

                X_a(ii) = sum( PhiywaveYptb_Row(:) * PhixRou_Row(:) )

            enddo
            enddo

            X_a = X_a + X_b

            !----------------------------------------------------------
            !   clear
            deallocate(Phiy,PhiyT)
            deallocate(Psi_A,Phiywave,PhiywaveT)
            deallocate(Phix,PhixT)
            deallocate(PhixRou_Row)
            deallocate(PhiywaveYptb_Row)
            return

        !------------------------------------------------------------------

    !!=======================================================================76
    end subroutine NLS5_4DVarCORE

    subroutine NLS5_EsbUpdate(      nx ,  ny ,  nz ,  nrens  ,  nobs  , DA_obs ,&
                            &       AA, Inflation, X_b ,  X_a)
    !!=======================================================================76
    !!
    !!  NLS5_EsbUpdate
    !!      samename:    AA(input)      ->      AA(output)
    !!
    !!=======================================================================76

        implicit none
        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   input
            integer,intent(in)          ::  nx
            integer,intent(in)          ::  ny
            integer,intent(in)          ::  nz
            integer,intent(in)          ::  nrens
            integer,intent(in)          ::  nobs
            type(DA_obs_bmc),intent(in) ::  DA_obs(nobs)
            real(8),intent(inout)       ::  AA(nx*ny*nz,nrens)
            real(8),intent(in)          ::  Inflation
            real(8),intent(in)          ::  X_b(nx*ny*nz)
            real(8),intent(in)          ::  X_a(nx*ny*nz)

            !   ptb
            real(8)                     ::  AAptb(nx*ny*nz,nrens)
            real(8)                     ::  ObsAAptb(nobs,nrens)
            
            !   parameter
            integer                     ::  ii, jj, kk
            integer                     ::  i,  j,  k

            !   singliar value decomposition
            real(8)                     ::  TT(nrens,nrens)
            real(8)                     ::  VV(nrens,nrens)
            real(8)                     ::  SigMa(nrens,nrens)
            real(8)                     ::  VVT(nrens,nrens)

            !   Psi_A parameter
            real(8)                     ::  Phiy(nobs,nrens)
            real(8)                     ::  PhiyT(nrens,nobs)
            real(8)                     ::  Psi_A(nrens,nrens)

            !   output
            real(8)                     ::  AAnew_ptb(nx*ny*nz,nrens)

            !----------------------------------------------------------

        !------------------------------------------------------------------

            !----------------------------------------------------------
            !   samples perturbation
            !           y' = Samples_y - Xb_y
            !           y'obs = Obs - ObsXb
            do ii=1,nobs
            do jj=1,nrens
                ObsAAptb(ii,jj) = DA_obs(ii)%obs_AA(jj) - DA_obs(ii)%obs_Xb
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
                PhiyT(jj,ii) = PhiyT(jj,ii) / (DA_obs(ii)%obs_err*DA_obs(ii)%obs_err)
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
    end subroutine NLS5_EsbUpdate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!===============================================================================84
End Module NLSCORE
