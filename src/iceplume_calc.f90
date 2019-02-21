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
! iceDepthK   - vertical grid index of grounding line depth         !
! plumeDepthK - vertical grid index of plume depth                  !
! osDepthK    - vertical grid index of overshoot depth              !
! areaP       - surface area of plume in contact with ice in        !
!               that cell [m^2]                                     !
! areaB       - surface area of plume in contact with ice out       !
!               that cell [m^2]                                     !
! meanVel     - ice tangental velocity [m s^-1]                     !
! depth       - calculated depth [m]                                !
! detr        - total detrainment volume [m^3 s^-1]                 !
! detrDz      - detrainment layer thickness [m]                     !
! mB          - background meltrate [m s^-1]                        !
!                                                                   !
! ==================================================================!
!
  integer, intent(in) :: ng, I
  real(r8), intent(in) :: dx, dy
  real(r8), intent(inout) :: fIni, tIni, sIni
  integer, intent(in) :: sgTyp
  real(r8), intent(in) :: sgDep, sgLen
  integer :: K, iTracer
  real(r8) :: RHO
!
  integer :: iceDepthK, plumeDepthK, osDepthK
  real(r8) :: areaP, areaB
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
  PLUME(ng) % rhoAm(I, K+1) = rhoAir
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
! Use entrainment plume model to calculate T, S, & flux in plume    !
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
        & +MAX(PLUME(ng) % mInt(I, K) - PLUME(ng) % mInt(I, K-1),       &
        &      0.0d0)
    ENDDO
!
! Separate entrainment / detrainment
!
    DO K = iceDepthK+1, N(ng)
      IF (PLUME(ng) % ent(I, K) .GT. 0.d0) THEN
        osDepthK = K
        detr = PLUME(ng) % ent(I, K)
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
! Use a detrainment model to distribute the detrainment plume in    !
! several layers. Details in iceplume_detrain.F.                    !
!                                                                   !
! *** This part is completely rewritten since Feb 20, 2019.         !
!                                                                   !
! ==================================================================!
!
  IF (fIni .GT. 0) THEN
!
! Find detrainment depth index, Calculate total detrainment volume
! and thickness
!
    plumeDepthK = -1
    detrDz = 1.d-7
    DO K = iceDepthK+1,osDepthK
      IF ((PLUME(ng) % rho(I, osDepthK) .LT.                            &
        &  PLUME(ng) % rhoAm(I, K)           ) .AND.                    &
        & (PLUME(ng) % rho(I, osDepthK) .GT.                            &
        &  PLUME(ng) % rhoAm(I, K+1)         )) THEN
        plumeDepthK = K
      ENDIF
    ENDDO
    PLUME(ng) % det(I, plumeDepthK) = detr
    PLUME(ng) % detI(I, plumeDepthK) = 1
    detrDz = detrDz + PLUME(ng) % dz(I, plumeDepthK)
!
! Use HALF parameterization to spread detrainment in several layers.
!
    CALL ICEPLUME_DETRAIN_HALF(ng, I,                                   &
                            &  iceDepthK, plumeDepthK, dx,              &
                            &  PLUME(ng) % lc(I, plumeDepthK),          &
                            &  detr, detrDz)
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
      & 0.5*(PLUME(ng) % s(I, osDepthK-1) +                             &
      &      PLUME(ng) % s(I, osDepthK))
    PLUME(ng) % trc(I, itemp) =                                         &
      & 0.5*(PLUME(ng) % t(I, osDepthK-1) +                             &
      &      PLUME(ng) % t(I, osDepthK))
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
      DO K = iceDepthK+1, osDepthK
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
        DO K = iceDepthK+1, osDepthK
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
