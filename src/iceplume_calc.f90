! ==================================================================!
!                                                                   !
! This subroutine call the iceplume core functions to calcualte     !
! plume status / melt rates / etc.                                  !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_CALC(ng, I, dx, dy, &
                       & QIni, TIni, SIni)
!
  USE mod_iceplume
  implicit none
!
! ==================================================================
! LOCAL VARIABLES:
! plumeAreaInCell :: surface area of plume in contact with ice in
! that cell (m^2)
! negSum, posSum :: sum of negative and positive contributions to
! the plume volume
! posNegRatio    :: ratio of the above
! meanVel :: ice tangental velocity
! ==================================================================
!
  integer, intent(in) :: ng, I
  real(r8), intent(in) :: dx, dy
  real(r8) :: QIni, TIni, SIni
  real(r8) :: rIni
  integer :: iceDepthK, plumeDepthK
  real(r8) :: negSum, posSum, posNegRatio
  real(r8) :: meanVel, depth
  real(r8) :: plumeAreaInCell
  real(r8) :: detr, detrDz
  integer :: K, iTracer
  real(r8) :: RHO
!
! Read in the subglacial discharge for this cell
!
  IF (useSheetPlume) THEN
      rIni = QIni/(wIni*dy)
  ELSEIF (useConePlume) THEN
      rIni = sqrt(2.d0*QIni/(pi*wIni))
  ELSE
      rIni = 0.d0
  ENDIF
!
! Create variables with temperature, salinity, density
! and velocity profiles for that column
!
  DO K = 1, N(ng)
      PLUME(ng) % dz(I, K) = &
        & PLUME(ng) % zW(I, K) - PLUME(ng) % zW(I, K-1)
  ENDDO
!
! If discharge input salinity is less or equal than zero, set it to
! a small number
!
  IF (SIni .LE. 0.d0) THEN
    SIni = 1.d-3
  ENDIF
!
! ==================================================================
! Find iceDepthK
! ==================================================================
!
  iceDepthK = -1
  IF (iceDepth .GE. 0) THEN
    iceDepthK = 0
  ELSE
    DO K = 0, N(ng)
      IF (PLUME(ng) % zW(I, K) .GE. iceDepth) THEN
        iceDepthK = K
        EXIT
      ENDIF
    ENDDO
  ENDIF
!
! ==================================================================
! Use plume model to calculate T, S, & flux in plume
! ==================================================================
!
  IF (QIni .GT. 0) THEN
!
! Run the plume model
!
    CALL ICEPLUME_PLUME_MODEL(ng, I, dy, iceDepthK, rIni, TIni, SIni)
!
! Calculate vertical plume volume flux
!
    DO k = 1, N(ng)
!
! After checking to see if we are above the base of the ice face...
!
      IF (K .GT. iceDepthK) THEN
!
! assuming specified plume horizontal extent (for sheet flow)...
!
        IF (useSheetPlume) THEN
          PLUME(ng) % volFlux(I, K) = &
            & (PLUME(ng) % r(I, K)) * (PLUME(ng) % w(I, K)) * dy
        ELSEIF (useConePlume) THEN
          PLUME(ng) % volFlux(I, K) = &
            & pi * ((PLUME(ng) % r(I, K))**2) * &
            & (PLUME(ng) % w(I, K)) / 2.
        ELSE
          PLUME(ng) % volFlux(I, K) = 0.0d0
        ENDIF
      ELSE
        PLUME(ng) % volFlux(I, K) = 0.0d0
      ENDIF
    ENDDO
!
! A couple of corrections:
! Even if plume is still buoyant, it cannot flow through the fjord surface
!
    PLUME(ng) % volFlux(I, N(ng)) = 0.0d0
!
! The initial volume flux is equal to runoff
!
    PLUME(ng) % volFlux(I, iceDepthK) = QIni
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
      PLUME(ng) % ent(I, K) = &
        & -(PLUME(ng) % volFlux(I, K) - PLUME(ng) % volFlux(I, K-1))
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
            IF (PLUME(ng) % ent(I, K) .GT. 0) &
              & PLUME(ng) % ent(I, K) = (PLUME(ng) % ent(I, K)) * posNegRatio
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
  ELSE  ! (QIni .EQ. 0)
!
! If no subglacial discharge, then there is no plume
!
    DO k = 1, N(ng)
      PLUME(ng) % r(I, K) = 0.d0
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
! ==================================================================
! Check if detrainment velocity meets the CFL criteria
!
! This section is added to deal with the surface intensification
! of ROMS grid. When the grid thickness is too thin, the output
! velocity is too high that violates the hroziontal advection CFL
! criteria. If you are confident that the detrainment velocity
! is small enough that would blow-up the model, turn checkCFL in
! mod_iceplume.F to .false.
! 2018/4/20 Chuning Wang
!
! Check if detrainment velocity meets Richardson number criteria
!
! This section is added to deal with the surface intensification
! of ROMS grid. When the grid thickness is too small, the output
! velocity is too high that violates the hroziontal advection CFL
! criteria. If you are confident that the detrainment velocity
! is small enough that would blow-up the model, turn checkCFL in
! mod_iceplume.F to .false.
! 2018/4/20 Chuning Wang
!
! To keep the code clean, this big chunk of code is wrapped in a
! subroutine ICEPLUME_CURVE_DET.
! ==================================================================
!
  IF (QIni .GT. 0) THEN
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
      CALL ICEPLUME_CURVE_DET(ng, I, &
                            & iceDepthK, plumeDepthK, &
                            & dx, dy, &
                            & detr, detrDz)
    ENDIF  ! checkCFL .OR. checkRiB
  ENDIF  ! Qini .GT. 0.0d0
!
! ==================================================================
! Calculate active and passive tracers
!
! Update log
! Add active tracers, T and S into this scheme. This is necessary
! when distributing detrainment in several layers (checkCFL), since
! rewriting T & S in PLUME directly is in general not a good idea.
! 2018/4/20 Chuning Wang
! ==================================================================
!
! Clear local plume tracer variables
!
  DO iTracer = 1, NT(ng)
    PLUME(ng) % trc(I, iTracer)    = 0.d0
    PLUME(ng) % trcCum(I, iTracer) = 0.d0
  ENDDO
!
  IF (QIni .GT. 0) THEN
!
! Write the active tracer (T & S) concentration to PLUME (ng) % trc
!
    PLUME(ng) % trc(I, isalt) = &
      & 0.5*(PLUME(ng) % s(I, plumeDepthK-1) + &
      &      PLUME(ng) % s(I, plumeDepthK))
    PLUME(ng) % trc(I, itemp) = &
      & 0.5*(PLUME(ng) % t(I, plumeDepthK-1) + &
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
          PLUME(ng) % trcCum(I, iTracer) = &
            & PLUME(ng) % trcCum(I, iTracer) + &
            & PLUME(ng) % trcIni(I, iTracer) * QIni
        ENDDO
      ENDIF
!
! Add up total sum of each tracer in plume
!
      DO K = iceDepthK+1, N(ng)
        IF (PLUME(ng) % ent(I, K) .LT. 0.) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcCum(I, iTracer) = &
              & PLUME(ng) % trcCum(I, iTracer) + &
              & (-PLUME(ng) % ent(I, K) * &
              &   PLUME(ng) % trcAm(I, K, iTracer))
          ENDDO
        ENDIF
      ENDDO
!
! Calculate concentration of tracer in outflow 
!
      DO iTracer = 3, NT(ng)
        PLUME(ng) % trc(I, iTracer) = &
          & PLUME(ng) % trcCum(I, iTracer) / detr
      ENDDO
    ENDIF
  ENDIF
!
! ==================================================================
! Calculate melt rates
! ==================================================================
!
  DO K = 1, N(ng)
!
! Check if we are above sea bed
!
    IF (K .LE. iceDepthK) THEN
!
! If not then there is no melting
!
      PLUME(ng) % mAv(I, K) = 0.d0
      PLUME(ng) % m(I, K) = 0.d0
      PLUME(ng) % mAm(I, K) = 0.d0
    ELSE
!
! ==================================================================
! If there is a plume in that cell, then need to calculate plume melt
! rate distinct to background melt rate. Plume melt rate is already
! encorporated in the plrume model, and taken into account in the
! temperature and salinity of the plume outflow. It is useful though
! to have it available as a diagnostic.
! ==================================================================
!
      plumeAreaInCell = 0.d0
      IF ((QIni .NE. 0) .AND. (useConePlume .OR. useSheetPlume)) THEN
        plumeAreaInCell = PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
!
        IF (plumeAreaInCell .GT. 0.0) THEN
          PLUME(ng) % m(I, K) = &
            & (PLUME(ng) % mInt(I, K) - &
            &  PLUME(ng) % mInt(I, K-1)) / &
            & plumeAreaInCell
        ELSE
          PLUME(ng) % m(I, K) = 0.d0
        ENDIF
      ELSE
!
! If no plume in that cell set plume melt rate to zero
!
        PLUME(ng) % m(I, K) = 0.d0
      ENDIF
!
! ==================================================================
! Calculate the background melt rate (i.e. not generated by plumes).
! This will then be used to update the temperature and salinity in the
! adjacent cells. Velocities are calculated at cell faces - find
! averages for cell centres. Does not include velocity perpendicular to
! ice - this differs depending on orientation of ice front
! ==================================================================
!
      IF (useBkgMelt) THEN
!
        meanVel = ((PLUME(ng) % vAm(I, K))**2 + &
                 & (PLUME(ng) % wAm(I, K))**2)**0.5d0
        depth = 0.5d0*(PLUME(ng) % zW(I, K-1) + PLUME(ng) % zW(I, K))
        CALL ICEPLUME_MELTRATE(PLUME(ng) % tAm(I, K), &
                             & PLUME(ng) % sAm(I, K), &
                             & meanVel, depth,     &
                             & PLUME(ng) % mAm(I, K), &
                             & PLUME(ng) % sB(I, K),  &
                             & PLUME(ng) % tB(I, K))
      ELSE
        PLUME(ng) % mAm(I, K) = 0.d0
      ENDIF
!
! Get average melt rate. This is useful for visualizing melt patterns
! and assessing overall melt rate of glacier.
! The following should apply to both conical and sheet plume models
!
      IF (QIni .NE. 0) THEN
        plumeAreaInCell = PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
        IF (plumeAreaInCell .LE. dy*PLUME(ng) % dz(I, K)) THEN
          IF (plumeAreaInCell .LE. 0) THEN
!
! If there is no plume in cell, then the melt rate is
! equal to the background melt rate.
!
            PLUME(ng) % mAv(I, K) = PLUME(ng) % m(I, K)
          ELSE
!
! If there is a plume in cell, calculate average melt rate
!
            PLUME(ng) % mAv(I, K) = &
              & (PLUME(ng) % m(I, K)*plumeAreaInCell + &
              &  PLUME(ng) % mAm(I, K) * &
              & (dy*(PLUME(ng) % dz(I, K))-plumeAreaInCell)) / &
              & (dy*(PLUME(ng) % dz(I, K)))
!
! Scale down background melt rate to account for area occupied by
! plume (necessary so that tendency terms aren't over estimated)
!
            PLUME(ng) % mAm(I, K) = &
              & PLUME(ng) % mAm(I, K) * &
              & (1 - plumeAreaInCell/(dy*(PLUME(ng) % dz(I, K))))
          ENDIF
        ELSE  ! plumeAreaInCell .GE. dy*dz(K)
!
! If the plume contact area is larger than the cell area, we
! assume there is no background melting
!
          PLUME(ng) % mAv(I, K) = &
            & (PLUME(ng) % m(I, K)) * plumeAreaInCell / &
            & (dy*(PLUME(ng) % dz(I, K)))
          PLUME(ng) % mAm(I, K) = 0.d0
        ENDIF
      ELSE
!
! If it is not a plume cell, then no plume melting.
!
        PLUME(ng) % m(I, K) = 0.d0
        PLUME(ng) % mAv(I, K) = PLUME(ng) % mAm(I, K)
      ENDIF  ! plume type
    ENDIF  ! above or below sea bed
  ENDDO
!
! ==================================================================
! Calculate thermodynamics
! ==================================================================
!
  IF (useBkgMelt) THEN
!
! Calculate the rate of temperature decrease in grid cells due to
! backgroud melting
!
    DO K = 1, N(ng)
!
! Check if above ice depth
!
      IF (K .LE. iceDepthK) THEN
!
! Below the ice depth, there is no heat and freshwater flux
!
        PLUME(ng) % fwFlux(I, K) = 0.d0
        PLUME(ng) % heatFlux(I, K) = 0.d0
        PLUME(ng) % tendT(I, K) = 0.d0
        PLUME(ng) % tendS(I, K) = 0.d0
      ELSE
!
! To convert from melt rate (m s^-1) to freshwater flux (kg m^-2 s^-1)
!
        PLUME(ng) % fwFlux(I, K) = &
          & -(PLUME(ng) % mAm(I, K)) * &
          & RHO(PLUME(ng) % tB(I, K), &
          &     PLUME(ng) % sB(I, K), &
          &     PLUME(ng) % zR(I, K))
!
! Latent heat required to melt that much ice (W m^-2)
!
        PLUME(ng) % heatFlux(I, K) = -(PLUME(ng) % fwFlux(I, K))*L
!
! If there is a plume, some of the freshwater is entrained into the
! plume. Scale fwFlux accordingly
!
!         plumeAreaInCell = &
!           & PLUME(ng) % a(I, K) - PLUME(ng) % a(I, K-1)
!         IF ((plumeAreaInCell .GT. 0) .AND. (correctMeltEnt)) THEN
!             PLUME(ng) % fwFlux(I, K) = &
!               & (PLUME(ng) % fwFlux(I, K))*(1-meltEnt)
!         ENDIF
!
! Compute tendencies (as for pkg/icefront in MITgcm)
!
        PLUME(ng) % tendT(I, K) = &
          & -(PLUME(ng) % heatFlux(I, K)) / c_i / rho_ref
        PLUME(ng) % tendS(I, K) = &
          & (PLUME(ng) % fwFlux(I, K)) * (PLUME(ng) % sAm(I, K)) / &
          & rho_ref
!
! Scale by icefrontlength, which is the ratio of the horizontal length
! of the ice front in each model grid cell divided by the grid cell
! area. (icefrontlength = dy / dxdy = 1 / dx)
!
        PLUME(ng) % tendT(I, K) = PLUME(ng) % tendT(I, K)/dx
        PLUME(ng) % tendS(I, K) = PLUME(ng) % tendS(I, K)/dx
      ENDIF
    ENDDO
  ELSE
    DO K = 1, N(ng)
      PLUME(ng) % fwFlux(I, K) = 0.d0
      PLUME(ng) % heatFlux(I, K) = 0.d0
      PLUME(ng) % tendT(I, K) = 0.d0
      PLUME(ng) % tendS(I, K) = 0.d0
    ENDDO
  ENDIF
END SUBROUTINE ICEPLUME_CALC
!
! ==================================================================
!
SUBROUTINE ICEPLUME_CURVE_DET(ng, I, &
                            & iceDepthK, plumeDepthK, &
                            & dx, dy, &
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
  real(r8) :: rho0, rhoUp, rhoDown, RiB, bvf
  real(r8) :: RHO
!
! For detrainment weight function
!
  real(r8) :: detrVel
  real(r8) :: detrWeight, detrWeightSum
  integer  :: KI, detrN, searchSwitch
  integer :: K, counter
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
        rhoUp = RHO(PLUME(ng) % tAm(I, K), &
                  & PLUME(ng) % sAm(I, K), &
                  & PLUME(ng) % zW(I, K))
      ELSE
        rhoUp = RHO(0.5*(PLUME(ng) % tAm(I, K+1) + &
                  &      PLUME(ng) % tAm(I, K)), &
                  & 0.5*(PLUME(ng) % sAm(I, K+1) + &
                  &      PLUME(ng) % sAm(I, K)), &
                  & PLUME(ng) % zW(I, K))
      ENDIF
!
      IF (K .EQ. 1) THEN
        rhoDown = RHO(PLUME(ng) % tAm(I, K), &
                    & PLUME(ng) % sAm(I, K), &
                    & PLUME(ng) % zW(I, K-1))
      ELSE
        rhoDown = RHO(0.5*(PLUME(ng) % tAm(I, K-1) + &
                    &      PLUME(ng) % tAm(I, K)), &
                    & 0.5*(PLUME(ng) % sAm(I, K-1) + &
                    &      PLUME(ng) % sAm(I, K)), &
                    & PLUME(ng) % zW(I, K-1))
      ENDIF
      rho0 = 0.5*(rhoUp + rhoDown)
      bvf = -g*(rhoUp-rhoDown)/(PLUME(ng) % dz(I, K)*rho0)
      minDetrDz2 = (detr**2*RiBmin/bvf/(dy**2))**0.25
    ELSE
      minDetrDz2 = 0.0d0
    ENDIF
!
    minDetrDz = MAX(minDetrDz, minDetrDz2)
!
! ==================================================================
! Distribute detrainmnet in several layers
!
! Update log
! Use a Gause function to smooth the distribution
! 2018/06/08 Chuning Wang
! ==================================================================
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
              IF ((PLUME(ng) % detI(I, K) .EQ. 0) .and. &
                & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
                KI = K
              ENDIF
            ENDDO
          ELSE
!
! If the plume has reached surface, search for one layer down
!
            DO K = 1,N(ng)-1
              IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .and. &
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
              IF ((PLUME(ng) % detI(I, K) .EQ. 0) .and. &
                & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
                KI = K
              ENDIF
            ENDDO
          ELSEIF (searchSwitch .EQ. -1) THEN
!
! Search one layer down
!
            DO K = 1,N(ng)-1
              IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .and. &
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
            detrWeight = &
              & EXP(-1.0 * &
              & (2 * (real(K)-real(plumeDepthK)) / real(detrN))**2)
            PLUME(ng) % det(I, K) = detrWeight * detrVel * dy * &
              & PLUME(ng) % dz(I, K)
          ENDIF
        ENDDO
!
! Normalize
!
        detrWeightSum = SUM(PLUME(ng) % det(I, :))
        DO K = 1,N(ng)
          PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K) * &
            & detr / detrWeightSum
        ENDDO
      ENDIF  ! detrN .GT. 1
END SUBROUTINE ICEPLUME_CURVE_DET
!
! ====== Thess subroutines are taken from MITgcm ===================
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
  rv =      a0 + (a1 + (a2 + a3*T)*T)*T &
  &     + (b0 + b1*T)*(S-sref) &
  &     + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
  &     + (  e0 + (e1 + e2*T)*T )*P*P
!
END
!
! =========================================================================
! Obsolete subroutines
! =========================================================================
!
! SUBROUTINE RHO_TO_W(profn, profr, profw)
! 
!     ! A quick program to extrapolate profiles from rho-points to w-points
!     implicit none
!     USE kinds
!     integer, intent(in) :: profn
!     real(r8) :: profr(profn), profw(profn+1)
! 
!     integer :: K
! 
!     profw(1) = profr(1)
!     profw(profn+1) = profr(profn)
!     DO K = 2, profn
!         profw(K) = 0.5d0*(profr(K-1)+profr(K))
!     ENDDO
! 
! END SUBROUTINE RHO_TO_W