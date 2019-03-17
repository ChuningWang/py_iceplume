!
! ==================================================================!
!                                                                   !
! These are the core functions of iceplume detrainment model.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN_FULL(ng, I,                                 &
                               & iceDepthK, plumeDepthK,                &
                               & dx, dy,                                &
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
! Local variables declaration
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
END SUBROUTINE ICEPLUME_DETRAIN_FULL
!
! ==================================================================!
!                                                                   !
! Use this function to calculate the thickness of detrainment.      !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN_HALF(ng, I,                                 &
                              &  iceDepthK, plumeDepthK,                &
                              &  dx, dy,                                &
                              &  detr, detrDz)
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
! Local variables declaration
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
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
! Check the horizontal advection CFL criteria to
! determine if to distribute the detrainment in several layers.
! Calculate maxinum velocity based on the CFL criteria
!
  maxVel = CuMax*dx/dt(ng)
  minDetrDz = detr/(maxVel*dy)
!
! Check the Richardson number criteria to determine
! if to distribute the detrainment in several layers.
!
! Calculate Brunt-Vaisala frequency
!
  K = plumeDepthK
  IF (K .EQ. N(ng)) THEN
    rhoUp = PLUME(ng) % rhoAm(I, K)
  ELSE
    rhoUp =0.5*(PLUME(ng) % rhoAm(I, K+1) +                             &
      &         PLUME(ng) % rhoAm(I, K))
  ENDIF
!
  IF (K .EQ. 1) THEN
    rhoDown = PLUME(ng) % rhoAm(I, K)
  ELSE
    rhoDown = 0.5*(PLUME(ng) % rhoAm(I, K-1) +                          &
      &            PLUME(ng) % rhoAm(I, K))
  ENDIF
  rho0 = PLUME(ng) % rhoAm(I, K)
  rhoP = 0.5*(PLUME(ng) % rho(I, K) + PLUME(ng) % rho(I, K-1))
  N2A = -g*(rhoUp-rhoDown)/(PLUME(ng) % dz(I, K)*rho0)
  N2P = -g*(rhoP -rho0   )/(PLUME(ng) % dz(I, K)*rho0)
  N2  = MAX(N2A, N2P, N2Bkg)
  minDetrDz2 = (detr**2*RiBmin/N2/dy**2)**0.25
!
  minDetrDz = MAX(minDetrDz, minDetrDz2)
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers                !
!                                                                   !
! Update log                                                        !
! Use a Gause function to smooth the distribution                   !
! 2018/06/08 Chuning Wang                                           !
!                                                                   !
! ==================================================================!
!
! Initiate searchSwitch for depthFinder == 2
!
  IF (depthFinder .EQ. 2) searchSwitch = 1
  counter = 0
  DO WHILE (detrDz .LT. minDetrDz)
    counter = counter + 1
!
! If detrDz is NaN, exit loop
!
    IF ((detrDz .EQ. (detrDz+1)) .OR. (counter .GT. 100)) EXIT
!
    IF (depthFinder .EQ. 1) THEN
!
! Search method 1, search one layer up until meet the surface
! Search for the nearest layer
!
      IF (PLUME(ng) % det(I, N(ng)) .EQ. 0.d0) THEN
!
! If the plume has not yet reached the surface, search for one layer
! up
!
        DO K = 2,N(ng)
          IF ((PLUME(ng) % detI(I, K) .EQ. 0) .and.                     &
            & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
            KI = K
          ENDIF
        ENDDO
      ELSE
!
! If the plume has reached surface, search for one layer down
!
        DO K = 1,N(ng)-1
          IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .and.                   &
            & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
            KI = K
          ENDIF
        ENDDO
      ENDIF
!
    ELSEIF (depthFinder .EQ. 2) THEN
!
! Search method 2, search one layer up then one layer down
! If the plume has reached surface, hard code searchSwitch to -1
!
      IF (PLUME(ng) % det(I, N(ng)) .GT. 0.d0) THEN
        searchSwitch = -1
      ENDIF
!
      IF (searchSwitch .EQ. 1) THEN
!
! Search one layer up
!
        DO K = 2,N(ng)
          IF ((PLUME(ng) % detI(I, K) .EQ. 0) .and.                     &
            & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
            KI = K
          ENDIF
        ENDDO
      ELSEIF (searchSwitch .EQ. -1) THEN
!
! Search one layer down
!
        DO K = 1,N(ng)-1
          IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .and.                   &
            & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
            KI = K
          ENDIF
        ENDDO
      ENDIF
      searchSwitch = searchSwitch*(-1)
!
    ELSE
!
! Other unknown search method will triger this error
!
      KI = -1
!
    ENDIF
!
! Update detrainment layer flag and detrainment depth
!
    PLUME(ng) % detI(I, KI) = 1
    detrDz = detrDz + PLUME(ng) % dz(I, KI)
!
! If KI is at the bottom of the plume, the search reaches the bottom and
! the barotropic velocity is possibly too large. Print the message to log
! file to alert the Users.
!
    IF (KI .EQ. iceDepthK) THEN
      write(*, *)  'ICEPLUME ALERT - KI = IceDepthK'
      write(*, *)  'checkCFL searched to the bottom of plume!'
      EXIT
    ENDIF
  ENDDO  ! WHILE detrDz .LT. minDetrDz
!
! Update detrainment volume flux
!
  detrN = SUM(PLUME(ng) % detI(I, :))
  IF (detrN .GT. 1) THEN
    detrVel = detr / (dy * detrDz)
    DO K = 1,N(ng)
      IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
        detrWeight = EXP(-1.0 * (2 *                                    &
          & (PLUME(ng) % zR(I, K)-                                      &
          &  PLUME(ng) % zR(I, plumeDepthK)) /                          &
          & detrDz)**2)
        PLUME(ng) % det(I, K) = detrWeight * detrVel * dy *             &
          & PLUME(ng) % dz(I, K)
      ENDIF
    ENDDO
!
! Normalize
!
    detrWeightSum = SUM(PLUME(ng) % det(I, :))
    DO K = 1,N(ng)
      PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K) *                   &
        & detr / detrWeightSum
    ENDDO
  ENDIF  ! detrN .GT. 1
END SUBROUTINE ICEPLUME_DETRAIN_HALF
