!
! ==================================================================!
!                                                                   !
! These are the core functions of  detrainment model.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN(ng, I,                                      &
                          & iceDepthK, plumeDepthK,                     &
                          & dx, dy,                                     &
                          & detr, detrDz)
  USE mod_iceplume
  implicit none
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK
  real(r8), intent(in) :: dx, dy
  real(r8), intent(in) :: detr
  real(r8), intent(inout) :: detrDz
!
! For checkCFL & checkRiB
!
  real(r8) :: maxVel, minDetrDz, minDetrDz2
  real(r8) :: rho0, rhoP, rhoUp, rhoDown, N2, N2A, N2P
!
! For detrainment weight function
!
  real(r8) :: detrVel
  real(r8) :: detrWeight, detrWeightSum
  integer  :: KI, detrN, searchSwitch
  integer :: K, counter
END SUBROUTINE ICEPLUME_DETRAIN
