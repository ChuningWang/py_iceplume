!
! ==================================================================!
!                                                                   !
! Use this function to calculate the velocity of detrainment.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN(ng, I,                                      &
                          & iceDepthK, plumeDepthK,                     &
                          & lc, detr)
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input/Output variables:                                           !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - glacier grounding line depth vertical grid index    !
! plumeDepthK - plume detrainment depth vertical grid index         !
! lc          - plume/ocean contact length [m]                      !
! detr        - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detrDz      - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
! maxVel      - maximum velocity based on Richardson number         !
!               criterion [m s^-1]                                  !
! Fr          - Froude number                                       !
! rhoP        - plume density [kg m^-3]                             !
! rho1        - upper layer density [kg m^-3]                       !
! rho2        - lower layer density [kg m^-3]                       !
! h1          - upper layer thickness [m]                           !
! h2          - lower layer thickness [m]                           !
! gRed        - reduced gravity [m s^-2]                            !
! isSurface   - logical, if the plume discharge at surface          !
!                                                                   !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
! potE        - potential energy of plume/ambient water             !
!                                                                   !
! detrVel     - detrainment velocity [m s^-1]                       !
! detrWeight  - detrainment weight function                         !
! detrSum     - detrainment volume flux before normalize [m^3 s^-1] !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK
  real(r8), intent(in) :: lc, detr
!
! Local variables declaration
!
  real(r8) :: detrDz, minDz, maxVel, Fr
  real(r8) :: rhoP, rho1, rho2, h1, h2, gRed, potE
  integer  :: udSwitch
  logical :: isSurface = .FALSE.
!
! For detrainment weight function
!
  real(r8) :: detrVel, detrWeight, detrSum
!
! Other local variables
!
  integer  :: KI, K, counter
  real(r8) :: rhoH
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
  detrDz = 0.0
!
! First, check if plume detrains in surface layer.
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
!
  IF (isSurface) THEN
!
! Use Froude number to limit detrainment flow speed.
!
    rhoP = PLUME(ng) % rho(I, plumeDepthK)
    rhoH = 0.0
    h1 = 0.0
    DO K = N(ng), 1, -1
      PLUME(ng) % detI(I, K) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, K)
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      rho1 = rhoH/h1
      gRed = MAX(g*(rho1-rhoP)/rhoRef, 0.0)
      Fr = detr/(lc*h1)/SQRT(gRed*h1)
      IF (Fr .LE. 1.0/RiB) EXIT
    ENDDO
  ELSE  ! (.NOT. isSurface)
!
! Use the Richardson number criteria to determine if to distribute
! the detrainment in several layers. First, calculate BV frequency.
!
    PLUME(ng) % detI(I, plumeDepthK) = 1
    detrDz = detrDz + PLUME(ng) % dz(I, plumeDepthK)
!
! Get density of the upper and lower layer.
!
    rhoH = 0.0
    h1 = 0.0
    DO K = plumeDepthK, N(ng)
      IF (K .EQ. plumeDepthK) THEN
        h1 = h1+0.5*PLUME(ng) % dz(I, K)
        rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h1 = h1+PLUME(ng) % dz(I, K)
        rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho1 = rhoH/h1
    rhoH = 0.0
    h2 = 0.0
    DO K = 1, plumeDepthK
      IF (K .EQ. plumeDepthK) THEN
        h2 = h2+0.5*PLUME(ng) % dz(I, K)
        rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h2 = h2+PLUME(ng) % dz(I, K)
        rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho2 = rhoH/h2
!
! Compute minimum detrainment thickness.
!
    gRed = MAX(-g*(rho1-rho2)/rhoRef, gRedBkg)
    maxVel = SQRT(gRed*h1*h2/(h1+h2))/RiB
    minDz = detr/maxVel/lc
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
    counter = 0
    KI = plumeDepthK
    potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
    DO WHILE ( (detrDz .LT. minDz) .AND. (counter .LT. 100) )
      counter = counter + 1
!
! If the plume has reached surface/bottom, set searchswitch to -1/1.
!
      IF     ( (PLUME(ng) % detI(I, N(ng)) .EQ. 1) .AND.                &
        &      (PLUME(ng) % detI(I, 0) .EQ. 1) ) THEN
        udSwitch = 0
        EXIT
      ELSEIF ( (PLUME(ng) % detI(I, N(ng)) .EQ. 1) .AND.                &
        &      (PLUME(ng) % detI(I, 0) .EQ. 0) ) THEN
        udSwitch = -1
      ELSEIF ( (PLUME(ng) % detI(I, N(ng)) .EQ. 0) .AND.                &
        &      (PLUME(ng) % detI(I, 0) .EQ. 1) ) THEN
        udSwitch = 1
      ENDIF
!
      IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
        DO K = 1, N(ng)
          IF ((PLUME(ng) % detI(I, K) .EQ. 0) .AND.                     &
            & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
            KI = K
          ENDIF
        ENDDO
      ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
        DO K = 1, N(ng)-1
          IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND.                   &
            & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
            KI = K
          ENDIF
        ENDDO
      ENDIF
!
! Search layer to minimize potential energy anomaly.
!
      potE = potE +                                                     &
        & g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
      IF (potE .GT. 0.0) THEN
        udSwitch = 1
      ELSE
        udSwitch = -1
      ENDIF
!
! Update detrainment layer flag and detrainment depth.
!
      PLUME(ng) % detI(I, KI) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, KI)
!
! This is another exit loop criteria.
!
      IF ( detrDz .EQ. (detrDz+1.0) ) EXIT
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers                !
!                                                                   !
! Update log                                                        !
! Use a Gause function to smooth the distribution.                  !
! 2018/06/08 Chuning Wang                                           !
!                                                                   !
! ==================================================================!
!
! Update detrainment volume flux
!
  detrVel = detr / (lc * detrDz)
  DO K = 1,N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % det(I, K) = detrVel * lc * PLUME(ng) % dz(I, K)
    ENDIF
  ENDDO
!
! Normalize
!
  detrSum = SUM(PLUME(ng) % det(I, :))
  DO K = 1,N(ng)
    PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K)*detr/detrSum
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN
