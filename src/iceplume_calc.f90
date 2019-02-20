!
! ==================================================================!
!                                                                   !
! This subroutine call the iceplume core functions to calcualte     !
! plume status / melt rates / etc.                                  !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_CALC(ng, I, dx, dy,                                 &
                       & fIni, tIni, sIni,                              &
                       & sgTyp, sgDep, sgLen)
!
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng      - grid identifier                                         !
! I       - point source index identifier                           !
! dx, dy  - grid length [m]                                         !
! fIni    - subglacial runoff volume flux [m^3 s^-1]                !
! tIni    - subglacial runoff temperature [degC]                    !
! sIni    - subglacial runoff salinity [PSU]                        !
! sgTyp   - subglacial runoff type identifier                       !
! sgDep   - subglacial runoff inject depth [m]                      !
! sgLen   - subglacial runoff discharge length [m]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! areaP            - surface area of plume in contact with ice in   !
!                   that cell [m^2]                                 !
! areaB            - surface area of plume in contact with ice out  !
!                   that cell [m^2]                                 !
! negSum, posSum   - sum of negative and positive contributions to  !
!                   the plume volume [m^3 s^-1]                     !
! posNegRatio      - ratio of the above                             !
! meanVel          - ice tangental velocity [m s^-1]                !
! depth            - calculated depth [m]                           !
! detr             - total detrainment volume [m^3 s^-1]            !
! detrDz           - detrainment layer thickness [m]                !
! mB               - background meltrate [m s^-1]                   !
!                                                                   !
! ==================================================================!
!
  integer, intent(in) :: ng, I
  real(r8), intent(in) :: dx, dy
  real(r8), intent(inout) :: fIni, tIni, sIni
  integer, intent(in) :: sgTyp
  real(r8), intent(in) :: sgDep, sgLen
  integer :: iceDepthK, plumeDepthK
  integer :: K, iTracer
  real(r8) :: RHO
!
  real(r8) :: areaP, areaB
  real(r8) :: negSum, posSum, posNegRatio
  real(r8) :: meanVel, depth
  real(r8) :: detr, detrDz
  real(r8) :: mB
!
! Calculate rho-layer depth, thickness, and ambient density
!
  DO K = 1, N(ng)
    PLUME(ng) % zR(I, K) =                                              &
      & 0.5d0 * (PLUME(ng) % zW(I, K-1) + PLUME(ng) % zW(I, K))
    PLUME(ng) % dz(I, K) =                                              &
      & PLUME(ng) % zW(I, K) - PLUME(ng) % zW(I, K-1)
    PLUME(ng) % rhoAm(I, K) =                                           &
      & RHO(PLUME(ng) % tAm(I, K),                                      &
      &     PLUME(ng) % sAm(I, K),                                      &
      &     PLUME(ng) % zR(I, K))
  ENDDO
!
! If discharge input salinity is less or equal than zero, set it to
! a small number
!
  IF (sIni .LE. 0.d0) THEN
    sIni = 1.d-3
  ENDIF
!
! ==================================================================!
!                                                                   !
! Find iceDepthK                                                    !
!                                                                   !
! ==================================================================!
!
  iceDepthK = -1
  IF (sgDep .GE. 0) THEN
    iceDepthK = 0
  ELSE
    DO K = 0, N(ng)
      IF (PLUME(ng) % zW(I, K) .GE. sgDep) THEN
        iceDepthK = K
        EXIT
      ENDIF
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! Use plume model to calculate T, S, & flux in plume                !
!                                                                   !
! ==================================================================!
!
  IF (fIni .GT. 0) THEN
!
! Run the plume model
!
    CALL ICEPLUME_ENTRAIN(ng, I, iceDepthK,                             &
                        & fIni, tIni, sIni,                             &
                        & sgTyp, sgLen)
!
! Calculate volume flux differential to give entrainment / extrainment
! First clear ent / det
!
    DO K = 1,N(ng)
      PLUME(ng) % ent(I, K) = 0.0d0
      PLUME(ng) % det(I, K) = 0.0d0
    ENDDO
!
    DO K = iceDepthK+1, N(ng)
      PLUME(ng) % ent(I, K) =                                           &
        & -(PLUME(ng) % f(I, K) - PLUME(ng) % f(I, K-1))                &
        & -MAX(PLUME(ng) % mInt(I, K) - PLUME(ng) % mInt(I, K-1),       &
        &      0.0d0)
    ENDDO
!
    IF (conserveMass) THEN
!
! Scale output to compensate for entrainment lost in expanding of
! output layer, i.e. so that there is no net flow over boundary
!
      negSum = 0.d0
      posSum = 0.d0
!
      DO K = iceDepthK+1, N(ng)
        IF ( PLUME(ng) % ent(I, K) .LT. 0 ) THEN
          negSum = PLUME(ng) % ent(I, K) + negSum
        ELSE
          posSum = PLUME(ng) % ent(I, K) + posSum
        ENDIF
      ENDDO
!
      IF ( negSum .NE. 0 ) THEN
        posNegRatio = -negSum / posSum
          DO K = 1,N(ng)
            IF (PLUME(ng) % ent(I, K) .GT. 0) THEN
              PLUME(ng) % ent(I, K) =                                   &
                & (PLUME(ng) % ent(I, K)) * posNegRatio
            ENDIF
        ENDDO
      ENDIF
    ENDIF
!
! Separate entrainment / detrainment
!
    DO K = iceDepthK+1, N(ng)
      IF (PLUME(ng) % ent(I, K) .GT. 0.d0) THEN
        PLUME(ng) % det(I, K) = PLUME(ng) % ent(I, K)
        PLUME(ng) % ent(I, K) = 0.d0
      ENDIF
    ENDDO
!
  ELSE  ! (fIni .EQ. 0)
!
! If no subglacial discharge, then there is no plume
!
    PLUME(ng) % f(I, 0) = 0.d0
    PLUME(ng) % w(I, 0) = 0.d0
    PLUME(ng) % t(I, 0) = 0.d0
    PLUME(ng) % s(I, 0) = 0.d0
    PLUME(ng) % a(I, 0) = 0.d0
    PLUME(ng) % mInt(I, 0) = 0.d0
    DO k = 1, N(ng)
      PLUME(ng) % f(I, K) = 0.d0
      PLUME(ng) % w(I, K) = 0.d0
      PLUME(ng) % t(I, K) = 0.d0
      PLUME(ng) % s(I, K) = 0.d0
      PLUME(ng) % a(I, K) = 0.d0
      PLUME(ng) % mInt(I, K) = 0.d0
      PLUME(ng) % ent(I, K) = 0.d0
      PLUME(ng) % det(I, K) = 0.d0
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! Check if detrainment velocity meets the CFL criteria              !
!                                                                   !
! This section is added to deal with the surface intensification    !
! of ROMS grid. When the grid thickness is too thin, the output     !
! velocity is too high that violates the hroziontal advection CFL   !
! criteria. If you are confident that the detrainment velocity      !
! is small enough that would blow-up the model, turn checkCFL in    !
! mod_iceplume.F to .false.                                         !
! 2018/4/20 Chuning Wang                                            !
!                                                                   !
! Check if detrainment velocity meets Richardson number criteria    !
!                                                                   !
! This section is added to deal with the surface intensification    !
! of ROMS grid. When the grid thickness is too small, the output    !
! velocity is too high that violates the hroziontal advection CFL   !
! criteria. If you are confident that the detrainment velocity      !
! is small enough that would blow-up the model, turn checkCFL in    !
! mod_iceplume.F to .false.                                         !
! 2018/4/20 Chuning Wang                                            !
!                                                                   !
! To keep the code clean, this big chunk of code is wrapped in a    !
! subroutine ICEPLUME_CURVE_DET.                                    !
!                                                                   !
! ==================================================================!
!
  IF (fIni .GT. 0) THEN
!
! Find detrainment depth index, Calculate total detrainment volume
! and thickness
!
    detr = 0.d0
    detrDz = 1.d-7
    DO K = 1,N(ng)
      IF (PLUME(ng) % det(I, K) .GT. 0.d0) THEN
        plumeDepthK = K
        detr = detr + PLUME(ng) % det(I, K)
        detrDz = detrDz + PLUME(ng) % dz(I, K)
        PLUME(ng) % detI(I, K) = 1
      ELSE
        PLUME(ng) % detI(I, K) = 0
      ENDIF
    ENDDO
!
    IF (checkCFL .OR. checkRiB) THEN
      CALL ICEPLUME_CURVE_DET(ng, I,                                    &
                            & iceDepthK, plumeDepthK,                   &
                            & dx, PLUME(ng) % lc(I, plumeDepthK),       &
                            & detr, detrDz)
    ENDIF  ! checkCFL .OR. checkRiB
  ELSE
    detr = 0.d0
    detrDz = 1.d-7
    plumeDepthK = iceDepthK+1
  ENDIF  ! fini .GT. 0.0d0
!
! ==================================================================!
!                                                                   !
! Calculate melt rates.                                             !
!                                                                   !
! ==================================================================!
!
! Initiate
!
  DO K = 1, N(ng)
    DO iTracer = 1, NT(ng)
      PLUME(ng) % trcB(I, K, iTracer)     = 0.d0
      PLUME(ng) % trcAmToB(I, K, iTracer) = 0.d0
    ENDDO
  ENDDO
!
  DO K = 1, N(ng)
!
! Check if we are above glacier grounding depth
!
    IF (K .LE. iceDepthK) THEN
!
! If not then there is no melting
!
      PLUME(ng) % m(I, K)  = 0.d0
      PLUME(ng) % mB(I, K) = 0.d0
!
    ELSE
!
! ==================================================================!
!                                                                   !
! If there is a plume in that cell, then need to calculate plume    !
! meltrate distinct to background melt rate. Plume melt rate is     !
! already encorporated in the plume model, and taken into account   !
! in the temperature and salinity of the plume outflow. It is       !
! useful though to have it available as a diagnostic.               !
!                                                                   !
! ==================================================================!
!
      areaP = PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
      areaP = MAX(areaP, 0.d0)
      areaB = PLUME(ng) % dz(I, K)*dy - areaP
      areaB = MAX(areaB, 0.d0)
!
      IF (areaP .GT. 0.0) THEN
        PLUME(ng) % m(I, K) = PLUME(ng) % mInt(I, K) -                  &
          &                   PLUME(ng) % mInt(I, K-1)
      ELSE
        PLUME(ng) % m(I, K) = 0.d0
      ENDIF
!
! ==================================================================!
!                                                                   !
! Calculate the background melt rate (not generated by plumes).     !
! This will then be used to update the temperature and salinity in  !
! the adjacent cells. Velocities are calculated at cell faces -     !
! find averages for cell centres. Does not include velocity         !
! perpendicular to ice - this differs depending on orientation of   !
! ice front.                                                        !
!                                                                   !
! ==================================================================!
!
      IF ((useBkgMelt) .AND. (areaB .GT. 0))THEN
!
        meanVel = SQRT((PLUME(ng) % vAm(I, K))**2 +                     &
          &            (PLUME(ng) % wAm(I, K))**2)
        CALL ICEPLUME_MELTRATE(PLUME(ng) % tAm(I, K),                   &
          &                    PLUME(ng) % sAm(I, K),                   &
          &                    meanVel,                                 &
          &                    PLUME(ng) % zR(I, K),                    &
          &                    mB,                                      &
          &                    PLUME(ng) % trcB(I, K, isalt),           &
          &                    PLUME(ng) % trcB(I, K, itemp))
!
! Calculate integrated background melt rate [m3 s-1]
!
        PLUME(ng) % mB(I, K) = mB*areaB
!
! Calculate equavilent temperature/salt/tracer turbulent flux
! [unit m^3 s^01] for background melting.
!
        PLUME(ng) % trcAmToB(I, K, itemp) =                             &
          & -PLUME(ng) % mB(I, K) *                                     &
          & (L + cI*(PLUME(ng) % trcB(I, K, itemp)-tIce))/cW
        PLUME(ng) % trcAmToB(I, K, isalt) =                             &
          & -PLUME(ng) % mB(I, K) *                                     &
          & (PLUME(ng) % trcB(I, K, isalt)-sIce)
        IF (useTracers) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcAmToB(I, K, iTracer) =                       &
              & -GamS * SQRT(CdBkg) * meanVel * areaB *                 &
              & PLUME(ng) % trcAm(I, K, iTracer)
          ENDDO
        ENDIF
!
      ELSE
        PLUME(ng) % mB(I, K) = 0.d0
      ENDIF
    ENDIF  ! above or below sea bed
  ENDDO
!
! ==================================================================!
!                                                                   !
! Calculate active and passive tracers                              !
!                                                                   !
! Update log                                                        !
! Add active tracers, T and S into this scheme. This is necessary   !
! when distributing detrainment in several layers (checkCFL), since !
! rewriting T & S in PLUME directly is in general not a good idea.  !
! 2018/4/20 Chuning Wang                                            !
!                                                                   !
! ==================================================================!
!
! Clear local plume tracer variables
!
  DO iTracer = 1, NT(ng)
    PLUME(ng) % trc(I, iTracer)    = 0.d0
    PLUME(ng) % trcCum(I, iTracer) = 0.d0
  ENDDO
!
  IF (fIni .GT. 0) THEN
!
! Write the active tracer (T & S) concentration to PLUME (ng) % trc
!
    PLUME(ng) % trc(I, isalt) =                                         &
      & 0.5*(PLUME(ng) % s(I, plumeDepthK-1) +                          &
      &      PLUME(ng) % s(I, plumeDepthK))
    PLUME(ng) % trc(I, itemp) =                                         &
      & 0.5*(PLUME(ng) % t(I, plumeDepthK-1) +                          &
      &      PLUME(ng) % t(I, plumeDepthK))
!
    IF (useTracers) THEN
!
! If use passive tracers, calculate passive tracer concentrations
! in the runoff
!
! Add ptracers in runoff
!
      IF (useInputTracers) THEN
        DO iTracer = 3, NT(ng)
          PLUME(ng) % trcCum(I, iTracer) =                              &
            & PLUME(ng) % trcCum(I, iTracer) +                          &
            & PLUME(ng) % trcIni(I, iTracer) * fIni
        ENDDO
      ENDIF
!
! Add up total sum of each tracer in plume
!
      DO K = iceDepthK+1, plumeDepthK
        IF (PLUME(ng) % ent(I, K) .LT. 0.) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcCum(I, iTracer) =                            &
              & PLUME(ng) % trcCum(I, iTracer) +                        &
              & (-PLUME(ng) % ent(I, K) *                               &
              &   PLUME(ng) % trcAm(I, K, iTracer))
          ENDDO
        ENDIF
      ENDDO
!
! Calculate concentration of tracer in outflow 
!
      DO iTracer = 3, NT(ng)
        PLUME(ng) % trc(I, iTracer) =                                   &
          & PLUME(ng) % trcCum(I, iTracer) / detr
      ENDDO
!
! If use melt water tracer, rewrite tracer concentration
!
      IF (useMeltTracers) THEN
!
! Calculate accumulated trcer amount from base to detrain depth
!
        PLUME(ng) % trcCum(I, NT(ng)-1) = 0.0d0
        PLUME(ng) % trcCum(I, NT(ng)  ) = 0.0d0
        DO K = iceDepthK+1, plumeDepthK
          PLUME(ng) % trcCum(I, NT(ng)-1) =                             &
            & PLUME(ng) % trcCum(I, NT(ng)-1) +                         &
            & (-PLUME(ng) % ent(I, K) *                                 &
            &   PLUME(ng) % trcAm(I, K, NT(ng)-1)) +                    &
            & PLUME(ng) % m(I, K) *                                     &
            & PLUME(ng) % trcIni(I, NT(ng)-1)
          PLUME(ng) % trcCum(I, NT(ng)  ) =                             &
            & PLUME(ng) % trcCum(I, NT(ng)  ) +                         &
            & (-PLUME(ng) % ent(I, K) *                                 &
            &   PLUME(ng) % trcAm(I, K, NT(ng))  )
        ENDDO
!
! Calculate tracer concentration
!
        PLUME(ng) % trc(I, NT(ng)-1) =                                  &
          & PLUME(ng) % trcCum(I, NT(ng)-1) / detr
        PLUME(ng) % trc(I, NT(ng)  ) =                                  &
          & PLUME(ng) % trcCum(I, NT(ng)  ) / detr
      ENDIF
    ENDIF
  ENDIF
!
! If use background melt tracer, rewrite tracer concentration
!
  IF ((useBkgMelt) .AND. (useTracers)) THEN
!
! Calculate tracer concentration in background melt boundary layer
!
    DO iTracer = 3, NT(ng)
      DO K = 1, N(ng)
        PLUME(ng) % trcB(I, K, iTracer) =                               &
          & -PLUME(ng) % trcAmToB(I, K, iTracer) /                      &
          & PLUME(ng) % mB(I, K)
      ENDDO
    ENDDO
    IF (useMeltTracers) THEN
      DO K = 1, N(ng)
        PLUME(ng) % trcB(I, K, NT(ng)) =                                &
          & PLUME(ng) % trcB(I, K, NT(ng)) +                            &
          & PLUME(ng) % trcIni(I, NT(ng))
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE ICEPLUME_CALC
!
! ==================================================================!
!                                                                   !
! Use this function to calculate melt rate.                         !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_MELTRATE(temp, salt, vel, depth, mdot, sB, tB)
!
  USE mod_iceplume
  implicit none
!
  real(r8) :: temp, salt, vel
  real(r8) :: depth, absVelocity
  real(r8) :: a, b, c
  real(r8), intent(inout) :: mdot, sB, tB
!
! Routine can't cope with zero velocity.
! Unlikely to occur anyway with currents, waves, convection etc.
! This isn't very physical, but will do for now.
!
  IF ( ABS(vel) .LT. velBkg ) vel = velBkg
!
! Calculate melt rate from 3 equation formualtion (as for plume
! models)
! Equations for Sb, Tb and mdot
!
  a = lambda1*(GamT*cW - GamS*cI)
  b = GamS*cI*(lambda1*salt - lambda2 - lambda3*depth +                 &
    &          tIce - (L/cI)) -                                         &
    &          GamT*cW*(temp - lambda2 - lambda3*depth + sIce)
  c = GamS*salt*(cI*(lambda2 + lambda3*depth - tIce) + L) +             &
    & GamT*sIce*cW*(temp - lambda2 - lambda3*depth)
!
  sB   = (1./(2.*a))*(-b - SQRT(b**2. - 4.*a*c))
  tB   = lambda1*sB + lambda2 + lambda3*depth
  mdot = GamS*SQRT(CdBkg)*ABS(vel)*(salt - sB)/sB
!
END SUBROUTINE ICEPLUME_MELTRATE
!
! ==================================================================!
!                                                                   !
! Use this function to calculate the thickness of detrainment.      !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_CURVE_DET(ng, I,                                    &
                            & iceDepthK, plumeDepthK,                   &
                            & dx, dy,                                   &
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
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
  IF (checkCFL) THEN
!
! If activiated, check the horizontal advection CFL criteria to
! determine if to distribute the detrainment in several layers.
!
! Calculate maxinum velocity based on the CFL criteria
!
    maxVel = CuMax*dx/dt(ng)
    minDetrDz = detr/(maxVel*dy)
  ELSE
    minDetrDz = 0.0d0
  ENDIF
!
  IF (checkRiB) THEN
!
! If activated, check the Richardson number criteria to determine
! if to distribute the detrainment in several layers.
!
! Calculate Brunt-Vaisala frequency
!
    K = plumeDepthK
    IF (K .EQ. N(ng)) THEN
      rhoUp = PLUME(ng) % rhoAm(I, K)
    ELSE
      rhoUp =0.5*(PLUME(ng) % rhoAm(I, K+1) +                           &
        &         PLUME(ng) % rhoAm(I, K))
    ENDIF
!
    IF (K .EQ. 1) THEN
      rhoDown = PLUME(ng) % rhoAm(I, K)
    ELSE
      rhoDown = 0.5*(PLUME(ng) % rhoAm(I, K-1) +                        &
        &            PLUME(ng) % rhoAm(I, K))
    ENDIF
    rho0 = PLUME(ng) % rhoAm(I, K)
    rhoP = 0.5*(PLUME(ng) % rho(I, K) + PLUME(ng) % rho(I, K-1))
    N2A = -g*(rhoUp-rhoDown)/(PLUME(ng) % dz(I, K)*rho0)
    N2P = -g*(rhoP -rho0   )/(PLUME(ng) % dz(I, K)*rho0)
    N2  = MAX(N2A, N2P, N2Bkg)
    minDetrDz2 = (detr**2*RiBmin/N2/dy**2)**0.25
  ELSE
    minDetrDz2 = 0.0d0
  ENDIF
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
END SUBROUTINE ICEPLUME_CURVE_DET
!
! ==================================================================!
!                                                                   !
! These subroutines are copied from MITgcm.                         !
!                                                                   !
! ==================================================================!
!
SUBROUTINE SW_TEMP(S, T, P, PR, rv)
!
! *=============================================================*
! | S/R  SW_TEMP
! | o compute in-situ temperature from potential temperature
! *=============================================================*
!
! REFERENCES:
! Fofonoff, P. and Millard, R.C. Jr
! Unesco 1983. Algorithms for computation of fundamental properties of
! seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
! Eqn.(31) p.39
! 
! Bryden, H. 1973.
! New Polynomials for thermal expansion, adiabatic temperature
! gradient and potential temperature of sea water.
! DEEP-SEA RES., 1973, Vol20,401-408.
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
! 
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S      :: salinity
! T      :: potential temperature
! P      :: pressure
! PR     :: reference pressure
! rv :: return value (in-situ temeparture in degree C)
!
  real(r8) :: S, T, P, PR
  real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
  CALL SW_PTMP  (S, T, PR, P, rv)
!
END
!
SUBROUTINE SW_PTMP  (S, T, P, PR, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_PTMP
! | o compute potential temperature as per UNESCO 1983 report.
! *=============================================================*
! \ev
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! ==================================================================
! SUBROUTINE SW_PTMP
!
! ==================================================================
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
! S  :: salinity    [psu      (PSS-78) ]
! T  :: temperature [degree C (IPTS-68)]
! P  :: pressure    [db]
! PR :: Reference pressure  [db]
! rv :: return value (potential temeparture in degree C)
!
  real(r8) :: S, T, P, PR
  real(r8) :: rv
!
! !LOCAL VARIABLES
! === local variables ===
!
  real(r8) :: del_P ,del_th, th, q
  real(r8) :: onehalf, two, three
  parameter ( onehalf = 0.5d0, two = 2.d0, three = 3.d0 )
  real(r8) :: adtg_val
!
! theta1
!
  del_P   = PR - P
  call sw_adtg(S, T, P, adtg_val)
  del_th  = del_P*adtg_val
  th      = T + onehalf*del_th
  q       = del_th
!
! theta2
!
  call sw_adtg(S, th, P+onehalf*del_P, adtg_val)
  del_th  = del_P*adtg_val
  th      = th + (1 - 1/sqrt(two))*(del_th - q)
  q       = (two-sqrt(two))*del_th + (-two+three/sqrt(two))*q
!
! theta3
!
  call sw_adtg(S, th, P+onehalf*del_P, adtg_val)
  del_th  = del_P*adtg_val
  th      = th + (1 + 1/sqrt(two))*(del_th - q)
  q       = (two + sqrt(two))*del_th + (-two-three/sqrt(two))*q
!
! theta4
!
  call sw_adtg(S, th, P+del_P, adtg_val)
  del_th  = del_P*adtg_val
  rv      = th + (del_th - two*q)/(two*three)
!
END
!
SUBROUTINE SW_ADTG  (S,T,P, rv)
!
! !DESCRIPTION: \bv
! *=============================================================*
! | S/R  SW_ADTG
! | o compute adiabatic temperature gradient as per UNESCO 1983 routines.
! *=============================================================*
! \ev
!
! started:
!          Armin Koehl akoehl@ucsd.edu
!
! !USES:
!
  USE mod_kinds
  IMPLICIT NONE
!
! === Global variables ===
!
! !INPUT/OUTPUT PARAMETERS:
! === Routine arguments ===
!
  real(r8) :: S,T,P
  real(r8) :: rv
!
! !LOCAL VARIABLES:
! === local variables ===
!
  real(r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  real(r8) :: sref
!
  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10
!
  b0 = +1.8932d-6
  b1 = -4.2393d-8
!
  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14
!
  d0 = -1.1351d-10
  d1 =  2.7759d-12
!
  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16
!
  rv =      a0 + (a1 + (a2 + a3*T)*T)*T                                 &
  &     + (b0 + b1*T)*(S-sref)                                          &
  &     + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P    &
  &     + (  e0 + (e1 + e2*T)*T )*P*P
!
END
