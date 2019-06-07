!
! ==================================================================!
!                                                                   !
! Use this function to calculate the velocity of detrainment.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN(ng, I,                                      &
     &                      iceDepthK, plumeDepthK, osDepthK,           &
     &                      lc, det)
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
! osDepthK    - plume overshoot depth vertical grid index           !
! lc          - plume/ocean contact length [m]                      !
! det         - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detDz       - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
! rhoP        - plume density [kg m^-3]                             !
! potE        - potential energy of plume/ambient water             !
!                                                                   !
! detVel      - detrainment velocity [m s^-1]                       !
! detWeight   - detrainment weight function                         !
! detSum      - detrainment volume flux before normalize [m^3 s^-1] !
!                                                                   !
! isSurface   - logical, if the plume discharge at surface          !
! isBottom    - logical, if the plume discharge at bottom           !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
!                                                                   !
! maxVel      - maximum velocity based on Richardson number         !
!               criterion [m s^-1]                                  !
! Fr          - Froude number                                       !
! rho1        - upper layer density [kg m^-3]                       !
! rho2        - lower layer density [kg m^-3]                       !
! h1          - upper layer thickness [m]                           !
! h2          - lower layer thickness [m]                           !
! gRed        - reduced gravity [m s^-2]                            !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, det
!
! Local variables declaration
!
  real(r8) :: detDz, minDz, detVel, detWeight, potE
  real(r8) :: rhoP
  integer  :: udSwitch
  logical  :: isSurface = .FALSE.
  logical  :: isBottom = .FALSE.
  real(r8) :: RHO
!
! For BV frequency formula
!
  real(r8) :: bvf0, bvf
!
! For detrainment weight function
!
  real(r8) :: detSum
!
! For calculation of internal wave speed
!
  real(r8) :: maxVel, Fr, Ri
  real(r8) :: rho1, rho2, h1, h2, gRed
  real(r8) :: rhoH, h01, h02
!
! Other local variables
!
  integer  :: KI, K, K1, K2, itrc, counter
  real(r8) :: cff, cff1, cff2, cff3, cff4
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
  IF (plumeDepthK .EQ. 1) isBottom = .TRUE.
  detDz = 0.0
  rhoP = RHO(PLUME(ng) % t(I, osDepthK), PLUME(ng) % s(I, osDepthK),    &
     &          PLUME(ng) % zR(I, plumeDepthK))
  IF (isSurface) THEN
    rhoH = 0.0
    h2 = 0.0
    DO K = iceDepthK+1, N(ng)
      h2 = h2+PLUME(ng) % dz(I, K)
      rhoH = rhoH + PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDDO
    rho2 = rhoH/h2
    gRed = MAX(-g*(rhoP-rho2)/rhoRef, 0.0)
    minDz = SQRT(det/(lc*gRed))
  ELSE
    KI = iceDepthK
    bvf0 = -g*(PLUME(ng) % rhoAm(I, KI+1) - PLUME(ng) % rhoAm(I, KI)) / &
        & (rhoRef*(PLUME(ng) % zW(I, KI+1) - PLUME(ng) % zW(I, KI)))
    DO K = iceDepthK+1, osDepthK-1
      bvf = -g*(PLUME(ng) % rhoAm(I, K+1) - PLUME(ng) % rhoAm(I, K)) /  &
        & (rhoRef*(PLUME(ng) % zW(I, K+1) - PLUME(ng) % zW(I, K)))
      IF ( bvf .GT. bvf0 ) THEN
        KI = K
        bvf0 = bvf
      ENDIF
    ENDDO
    rhoH = 0.0
    h1 = 0.0
    DO K = KI+1, N(ng)
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH + PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDDO
    rho1 = rhoH/h1
!
    rhoH = 0.0
    h2 = 0.0
    DO K = 1, KI 
      h2 = h2+PLUME(ng) % dz(I, K)
      rhoH = rhoH + PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDDO
    rho2 = rhoH/h2
!
    gRed = -g*(rho1-rho2)/rhoRef
    KI = MAXLOC(PLUME(ng) % w(I, :), 1)
    cff1 = PLUME(ng) % f(I, KI)/PLUME(ng) % w(I, KI)/lc
    Ri = gRed*cff1/(PLUME(ng) % w(I, KI)**2)
    IF (Ri .LT. 6.0) THEN
      cff = (0.7*Ri**0.17)*PLUME(ng) % w(I, KI)
    ELSE
      cff = 0.95*PLUME(ng) % w(I, KI)
    ENDIF
    minDz = det/lc/cff
  ENDIF
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
  counter = 0
  KI = plumeDepthK
  PLUME(ng) % detI(I, KI) = 1
  PLUME(ng) % detFrac(I, KI) = 1.0
  detDz = PLUME(ng) % dz(I, KI)
  K1 = KI
  K2 = KI
  potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
  IF (potE .GT. 0.0) THEN
    udSwitch = 1
  ELSE
    udSwitch = -1
  ENDIF
  DO WHILE ( (detDz .LT. minDz) .AND. (counter .LT. 100) )
    counter = counter + 1
    IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
      DO K = 1, N(ng)
        IF ((PLUME(ng) % detI(I, K  ) .EQ. 0) .AND.                     &
     &      (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
          KI = K
        ENDIF
      ENDDO
      K1 = KI
    ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
      DO K = 1, N(ng)-1
        IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND.                     &
     &      (PLUME(ng) % detI(I, K  ) .EQ. 0)) THEN
          KI = K
        ENDIF
      ENDDO
      K2 = KI
    ENDIF
    IF (KI .EQ. N(ng)) isSurface = .TRUE.
    IF (KI .EQ. 1) isBottom = .TRUE.
!
! Update detrainment layer flag and detrainment depth.
!
    PLUME(ng) % detI(I, KI) = 1
    detDz = detDz + PLUME(ng) % dz(I, KI)
    PLUME(ng) % detFrac(I, KI) = 1.0
!
! Determine next layer to minimize potential energy anomaly.
!
    potE = potE +                                                       &
     &     g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
!
! IF layer reached surface/bottom, force udSwitch to -1/1
!
    IF     ( (       isSurface ) .AND. ( .NOT. isBottom ) ) THEN
      udSwitch = -1
    ELSEIF ( ( .NOT. isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 1
    ELSEIF ( (       isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 0
    ENDIF
!
! Other exit loop criteria.
!
    IF ( (detDz .EQ. (detDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
  ENDDO
!
! Only take partial grid from top/bottom layer to satisfy potE=0
!
  IF ( (.NOT. isSurface) .AND. (counter .GT. 0) ) THEN
    IF ( potE .GT. 0.0 ) THEN
      KI = K2
    ELSE
      KI = K1
    ENDIF
    cff = 1 - potE / (g*(PLUME(ng) % rhoAm(I, KI)-rhoP)) /              &
     &        PLUME(ng) % dz(I, KI)
    IF ( (cff .GE. 0) .AND. (cff .LT. 1) ) THEN
      PLUME(ng) % detFrac(I, KI) = cff
      detDz = detDz - PLUME(ng) % dz(I, KI)*(1 - cff)
    ENDIF
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
  detVel = det / (lc * detDz)
  DO K = 1,N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
!       detWeight = EXP(-1.0 *                                          &
!      &    ((PLUME(ng) % zR(I, K)-                                     &
!      &      PLUME(ng) % zR(I, plumeDepthK)) /                         &
!      &     detDz)**2)
      detWeight = 1.0
!
! Distribute the detrainment.
!
      PLUME(ng) % det(I, K) = detWeight * detVel * lc *                 &
     &  PLUME(ng) % detFrac(I, K) * PLUME(ng) % dz(I, K)
    ENDIF
  ENDDO
!
! Normalize
!
  detSum = SUM(PLUME(ng) % det(I, :))
  DO K = 1,N(ng)
    PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K)*det/detSum
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN
