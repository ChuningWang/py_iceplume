! ==================================================================!
!                                                                   !
! These are the core functions of ICEPLUME model.                   !
!                                                                   !
! ==================================================================!
!
! ==================================================================
! Program to calculate the plume shape, S, T, V etc.
! ==================================================================
!
SUBROUTINE iceplume_plume_model(ng, rIni, tIni, sIni)
!
  USE mod_iceplume
  implicit none
  integer, intent(in) :: ng
  real(r8), intent(in) :: rIni, tIni, sIni
!
! ==================================================================
! Local variables for ODEPACK
! ==================================================================
!
  integer :: IOPT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW, MF, NEQ
  real(r8), parameter :: RTOL = 1.0d-5
  real(r8), parameter :: ATOL = 1.0d-5
  real(r8) :: RWORK(116), Y(6)
  real(r8) :: zIn, zOut
!
! Y is input/output vector for DLSODE
!   Y(1) = plume thickness/radius
!   Y(2) = plume velocity
!   Y(3) = plume temperature
!   Y(4) = plume salinity
!   Y(5) = plume area
!   Y(6) = area integrated melt
!
! ==================================================================
! Other local variables
! ==================================================================
!
  real(r8) :: RHO, temperature, salinity, depth
  real(r8) :: tAmbient, sAmbient
  real(r8) :: rhoPlume, rhoAmbient
  real(r8) :: vAmbient, wAmbient, meanVel, meltAmbient
  real(r8) :: plumeAreaInCell, plumeMass, meltMass
  real(r8) :: sB, tB
!
! Plume models
!
  external JENKINS, HALFCONE, JEX
!
! ==================================================================
! For ODEPACK solver. See ODEPACK documentation and source code in
! Cowton et al. 2015.
! ==================================================================
!
  NEQ = 6
  LRW = 116
  LIW = 116
!
  ITOL = 1
  ITASK = 1
  ISTATE = 1
  IOPT = 0
  MF = 10
  IWORK(7) = 2  ! To limit number of times repeat error messages are printed
!
! ==================================================================
! Initial conditions
! ==================================================================
!
  Y(1) = rIni  ! initial plume thickness (sheet model)
               ! or radius (halfcone model)
  Y(2) = wIni  ! initial vertical velocity
  Y(3) = tIni  ! initial temperature
  Y(4) = sIni  ! initial salinity
  Y(5) = 0.0   ! integrated contact area
  Y(6) = 0.0   ! integrated melt rate
!
! Prepare profiles
!
  DO K = 0, Nr
    PLUME(ng) % r(K) = 0.0
    PLUME(ng) % w(K) = 0.0
    PLUME(ng) % t(K) = 0.0
    PLUME(ng) % s(K) = 0.0
    PLUME(ng) % a(K) = 0.0
    PLUME(ng) % mInt(K) = 0.0
    PLUME(ng) % rho(K) = 0.0
  ENDDO
  DO K = 1, Nr
    PLUME(ng) % zR(K) = 0.5d0 * (PLUME(ng) % zW(K-1) + &
      &                          PLUME(ng) % zW(K))
  ENDDO
!
  plumeDepthK = 0
  rhoPlume = 0.d0
!
! Start at bottom of ice face
!
  zIn = PLUME(ng) % zW(iceDepthK)
!
! Next point at which to retrieve values
!
  zOut = PLUME(ng) % zW(iceDepthK+1)
!
! Set initial conditions
!
  PLUME(ng) % r(iceDepthK) = Y(1)
  PLUME(ng) % w(iceDepthK) = Y(2)
  PLUME(ng) % t(iceDepthK) = Y(3)
  PLUME(ng) % s(iceDepthK) = Y(4)
  PLUME(ng) % a(iceDepthK) = Y(5)
  PLUME(ng) % mInt(iceDepthK) = Y(6)
!
! ==================================================================
! Move up through water column from lowest layer
! ==================================================================
!
  DO K = iceDepthK+1, Nr
!
! Replicate current depth
!
  depth = zIn
!
! ==================================================================
! Use DLSODE to solve plume properties.
! ==================================================================
!
! Check to make sure plume hasn't reached neutral buoyancy in a lower
! layer
!
    IF (ISTATE .GT. -1) THEN
      IF (useSheetPlume) THEN
        CALL DLSODE (JENKINS, NEQ, Y, zIn,   &
          & zOut, ITOL, RTOL, ATOL, ITASK,   &
          & ISTATE, IOPT, RWORK, LRW, IWORK, &
          & LIW, JEX, MF)
      ELSEIF (useConePlume) THEN
        CALL DLSODE (HALFCONE, NEQ, Y, zIn,  &
          & zOut, ITOL, RTOL, ATOL, ITASK,   &
          & ISTATE, IOPT, RWORK, LRW, IWORK, &
          & LIW, JEX, MF)
      ENDIF
!
! ==================================================================
! Test to see if neutral buoyancy has now been reached. If solver
! returns ISTATE = -1, then it has been unable to meet required
! tolerances at this level. This generally occurs because plume
! has reached neutral buoyancy and run out of momentum, and so is no
! longer rising. At this point, we therefore end the call to the
! plume model. Our aim is to catch the plume at the point of neutral
! buoyancy. We therefore perform a manual comparrison of ambient and
! plume density. If plume density >= ambient density we assign
! ISTATE = -1, again ending the call to the plume model.
! ==================================================================
!
! Calculate plume density (rho = RHO(temp, salt, depth))
!
      rhoPlume = RHO(Y(3), Y(4), depth)
!
! Calculate ambient density
!
      IF (K .EQ. 0) THEN
        tAmbient = PLUME(ng) % tAm(1)
        sAmbient = PLUME(ng) % sAm(1)
      ELSEIF (K .EQ. Nr) THEN
        tAmbient = PLUME(ng) % tAm(Nr)
        sAmbient = PLUME(ng) % sAm(Nr)
      ELSE
        tAmbient = .5*(PLUME(ng) % tAm(K) + PLUME(ng) % tAm(K+1))
        sAmbient = .5*(PLUME(ng) % sAm(K) + PLUME(ng) % sAm(K+1))
      ENDIF
      rhoAmbient = RHO(tAmbient, sAmbient, depth)
      write(*, *)  K, depth, rhoAmbient, rhoPlume
!
      IF ((rhoPlume .GT. rhoAmbient) .OR. &
        & (K .EQ. Nr)) THEN
        ISTATE = -1
      ENDIF
!
! If ISTATE is now < 0, then plume has reached neutral buoyancy 
!
      IF (ISTATE .LT. 0) THEN
!
! If we have reached neutral buoyancy then there is no volume flux out
! of this cell, so plume area and velocity equal zero. Other values are
! kept for use in determining plume outflow properties.
!
        Y(1) = 0.d0
        Y(2) = 0.d0
      ELSE
!
! If the plume has not reached neutral buoyancy, then we assign a depth
! at which to calculate the next value and loop round to call the plume
! model again. Make sure we're not at the surface
!
        IF (K .NE. Nr) THEN
          zIn = zOut
          zOut = PLUME(ng) % zW(K+1)
        ENDIF
      ENDIF
!
! ==================================================================
! Make corrections for background melt water entrainment.
! (halfcone model only)
! ==================================================================
!
! Assuming some portion of the melt water is entrained into the plume,
! the acutual temperature and salinity should be lower in plume due to
! this extra entrainment. This section calculates background melting,
! calculates melt water entrainment and updates T & S in plume using a
! conservative mixing scheme.
!
      IF ((useConePlume) .AND. (correctMeltEnt)) THEN
!
! Check if plume has risen / extends to the whole grid
!
        IF ((K .GT. 0) .AND. (2.0*Y(1) .LT. dy)) THEN
          tAmbient = PLUME(ng) % tAm(K)
          sAmbient = PLUME(ng) % sAm(K)
          vAmbient = PLUME(ng) % vAm(K)
          wAmbient = PLUME(ng) % wAm(K)
          meanVel = (vAmbient**2 + wAmbient**2)**0.5
!
! Calculate background melt rate [m s^-1]
!
          CALL iceplume_meltrate(tAmbient,       &
                               & sAmbient,       &
                               & meanVel, depth, &
                               & meltAmbient, sB, tB)
!
! Calculate freshwater flux from melt water [kg s^-1]
!
          plumeAreaInCell = Y(5) - PLUME(ng) % a(K-1)
          meltMass = &
            & meltEnt * meltAmbient * rho_ice * &
            & (dy*(PLUME(ng) % dz(K)) - plumeAreaInCell)
!
! Calculate plume mass [kg s^-1]
!
          plumeMass = 0.5 * pi * Y(1)**2 * Y(2) * rhoPlume
!
! Update plume temperature and salinity with conservative mixing
!
          Y(3) = &
            & (Y(3)*plumeMass+tB*meltMass) / &
            & (plumeMass+meltMass)
          Y(4) = &
            & (Y(4)*plumeMass+sB*meltMass) / &
            & (plumeMass+meltMass)
        ENDIF
      ENDIF
!
    ELSE  ! (ISTATE .LE. -1)
!
! This section is entered once the plume has reached neutral buoyancy
! once plume has reached neutral buoyancy, no plume values
!
      Y(1) = 0.0
      Y(2) = 0.0
      Y(3) = 0.0
      Y(4) = 0.0
      Y(5) = 0.0
      Y(6) = 0.0
    ENDIF
!
! ==================================================================
! Save results
! ==================================================================
!
    PLUME(ng) % r(K) = Y(1)
    PLUME(ng) % w(K) = Y(2)
    PLUME(ng) % t(K) = Y(3)
    PLUME(ng) % s(K) = Y(4)
    PLUME(ng) % a(K) = Y(5)
    PLUME(ng) % mInt(K) = Y(6)
!
!    write(*, *)  K, &
!      & PLUME(ng) % zW(K), &
!      & PLUME(ng) % s(K), sAmbient, &
!      & PLUME(ng) % t(K), tAmbient, &
!      & PLUME(ng) % r(K)
!
  ENDDO
END
!
! ==================================================================!
!                                                                   !
! Use this function to calculate melt rate.                         !
!                                                                   !
! ==================================================================!
!
SUBROUTINE iceplume_meltrate(temperature, salinity, &
                           & velocity, depth, &
                           & mdot, Sb, Tb)
!
  USE mod_iceplume
  implicit none
!
  real(r8) :: temperature, salinity, velocity
  real(r8) :: depth, absVelocity
  real(r8) :: a, b, c
  real(r8), intent(inout) :: mdot, Sb, Tb
!
! Routine can't cope with zero velocity.
! Unlikely to occur anyway with currents, waves, convection etc.
! This isn't very physical, but will do for now.
!
  IF ( velocity .LT. backgroundVel ) velocity = backgroundVel
!
  absVelocity = abs(velocity)
!
! Calculate melt rate from 3 equation formualtion (as for plume models)
! Equations for Sb, Tb and mdot
!
  a = lambda1*(GamT*c_w-GamS*c_i)
!
  b = GamS*c_i*(lambda1*salinity-lambda2-lambda3*depth+ &
   &         iceTemp-(L/c_i)) &
   &        -GamT*c_w*(temperature-lambda2-lambda3*depth)
!
  c = GamS*salinity*(c_i*(lambda2+lambda3*depth-iceTemp)+L)
!
  Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5))      ! Sb
  Tb   = lambda1*Sb+lambda2+lambda3*depth            ! Tb
  mdot = GamS*(Cd**0.5)*absVelocity*(salinity-Sb)/Sb ! mdot
!
END
!
! ==================================================================!
!                                                                   !
! Plume models. These functions come from Dr. Tom Cowton.           !
!                                                                   !
! ==================================================================!
!
SUBROUTINE  HALFCONE (NEQ, T, Y, YDOT)
!
  USE mod_iceplume
  implicit none
!
  integer :: NEQ
  real(r8) :: T, Y(6), YDOT(6)
  real(r8) :: Tambient, Sambient, rho_0, rho_1
  real(r8) :: mdot, Sb, Tb
  real(r8) :: a, b, c
  real(r8) :: RHO
!
! Interpolate from imposed ambient profiles
!
  IF (T .LE. PLUME(ngr) % zR(1)) THEN
    Tambient = PLUME(ngr) % tAm(1)
    Sambient = PLUME(ngr) % sAm(1)
  ELSEIF (T .GE. PLUME(ngr) % zR(Nr)) THEN
    Tambient = PLUME(ngr) % tAm(Nr)
    Sambient = PLUME(ngr) % sAm(Nr)
  ELSE
    CALL LININT(Nr, PLUME(ngr) % zR, PLUME(ngr) % tAm, T, Tambient)
    CALL LININT(Nr, PLUME(ngr) % zR, PLUME(ngr) % sAm, T, Sambient)
  ENDIF
!
  rho_1 = RHO(Y(3), Y(4), T)
  rho_0 = RHO(Tambient, Sambient, T)
!
! Equations for Sb, Tb and mdot
!
  a = lambda1*(GamT*c_w-GamS*c_i)
!
  b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*T +  &
      &         iceTemp-(L/c_i)) -                &
      &         GamT*c_w*(Y(3)-lambda2-lambda3*T)
!
  c = GamS*Y(4)*(c_i*(lambda2+lambda3*T-iceTemp)+L)
!
  Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) ! Sb
  Tb   = lambda1*Sb+lambda2+lambda3*T           ! Tb
  mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb       ! mdot
!
! Differential equations
! Plume radius
!
  YDOT(1) = 2.*E_0+4.*mdot/(pi*Y(2)) - &
          & Y(1)*g*(rho_0-rho_1)/(2.*Y(2)*Y(2)*rho_ref)+2.*Cd/pi
!
! Plume vertical velocity
!
  YDOT(2) = -2.*E_0*Y(2)/Y(1)-4.*mdot/(pi*Y(1))+g* &
      &          (rho_0-rho_1)/(Y(2)*rho_ref)-4.*Cd*Y(2)/(pi*Y(1))
!
! Plume temperature
!
  YDOT(3) = 2.*E_0*(TAMBIENT-Y(3))/Y(1)+4.*mdot* &
      &           (Tb-Y(3))/(pi*Y(1)*Y(2))-4.* &
      &           GamT*(Cd**0.5)*(Y(3)-Tb)/(pi*Y(1))
!
! Plume salinity
!
  YDOT(4) = 2.*E_0*(Sambient-Y(4))/Y(1)+4.*mdot*(Sb-Y(4))/ &
      &          (pi*Y(1)*Y(2))-4.*GamS*(Cd**0.5)*(Y(4)-Sb)/(pi*Y(1))
!
! Along-plume integrated contact area and melt rate
!
  YDOT(5) = 2.*Y(1)
  YDOT(6) = 2.*Y(1)*mdot
!
END SUBROUTINE HALFCONE
!
! =========================================================================
!
SUBROUTINE  JENKINS (NEQ, T, Y, YDOT)
!
  USE mod_iceplume
  implicit none
!
  integer ::  NEQ
  real(r8) :: T, Y(6), YDOT(6)
  real(r8) :: Tambient, Sambient, rho_0, rho_1
  real(r8) :: mdot, Sb, Tb
  real(r8) :: a,b,c
  real(r8) :: RHO
!
! Interpolate from imposed ambient profiles
!
  IF (T .LE. PLUME(ngr) % zR(1)) THEN
    Tambient = PLUME(ngr) % tAm(1)
    Sambient = PLUME(ngr) % sAm(1)
  ELSEIF (T .GE. PLUME(ngr) % zR(Nr)) THEN
    Tambient = PLUME(ngr) % tAm(Nr)
    Sambient = PLUME(ngr) % sAm(Nr)
  ELSE
    CALL LININT(Nr, PLUME(ngr) % zR, PLUME(ngr) % tAm, T, Tambient)
    CALL LININT(Nr, PLUME(ngr) % zR, PLUME(ngr) % sAm, T, Sambient)
  ENDIF
!
  rho_1 = RHO(Y(3), Y(4), T)
  rho_0 = RHO(Tambient, Sambient, T)
!
! Equations for Sb, Tb and mdot
!
  a = lambda1*(GamT*c_w-GamS*c_i)
!
  b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*T+ &
      &         iceTemp-(L/c_i)) &
      &        -GamT*c_w*(Y(3)-lambda2-lambda3*T)
!
  c = GamS*Y(4)*(c_i*(lambda2+lambda3*T-iceTemp)+L)
!
  Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) ! Sb
  Tb   = lambda1*Sb+lambda2+lambda3*T           ! Tb
  mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb       ! mdot
!
! Differential equations
! Plume thickness
!
  YDOT(1)=2*E_0+Cd-(g*Y(1)/(Y(2)**2))*(rho_0-rho_1) &
  &  /rho_ref+2*mdot/Y(2)
!
! Plume vertical velocity
!
  YDOT(2)=-(Y(2)/Y(1))*(E_0+Cd+mdot/Y(2)) &
  &  +(g/Y(2))*(rho_0-rho_1)/rho_ref
!
! Plume temperature
!
  YDOT(3)=E_0*Tambient/Y(1)-(Y(3)/Y(1)) &
  &  *(E_0+mdot/Y(2))+(mdot/(Y(1)*Y(2))) &
  &  *(Tb-(L/c_w)-(c_i/c_w)*(Tb-iceTemp))
!
! Plume salinity
!
  YDOT(4)=E_0*Sambient/Y(1)-(Y(4)/Y(1)) &
  &  *(E_0+mdot/Y(2))
!
! along-plume integrated contact area and melt rate
!
  YDOT(5) = dy  ! This is constant in sheet model
  YDOT(6) = dy * mdot
!
END
!
! =========================================================================
!
DOUBLE PRECISION FUNCTION RHO(T,S,z)
!
! Equation of state (UNESCO 1983)
!     T = temperature (deg C)
!     S = salinity (PSU)
!     z = depth (m)
!
  DOUBLE PRECISION T,S,z
  DOUBLE PRECISION rho_0, g, P
  DOUBLE PRECISION kw, Aw, Bw, k0
  DOUBLE PRECISION bulk_modulus
  DOUBLE PRECISION A, B, C, rho_w,rho_zero
!
  PARAMETER(rho_0=1027)
  PARAMETER(g=9.81)
!
  P= rho_0*g*abs(z)*1.0E-5
!
! RHO_1 (in situ)
  kw= 19652.21+ 148.4206*T- 2.327105*T**2+ &
   &    1.360477e-2*(T**3)-5.155288e-5*(T**4)
  Aw= 3.239908+ 1.43713e-3*T+ 1.16092e-4*T**2- &
   &    5.77905e-7*T**3
  Bw= 8.50935e-5- 6.12293e-6*T + 5.2787e-8*(T**2)
  k0= kw + (54.6746- 0.603459*T+ 1.09987e-2*(T**2) &
   &    -6.1670e-5*(T**3))*S +(7.944e-2 + 1.6483e-2* &
   &    T- 5.3009e-4*(T**2))*(S**1.5)
  A=  Aw+ (2.2838e-3- 1.0981e-5*T- 1.6078e-6*(T**2)) &
   &    *S+ 1.91075e-4*(S**1.5)
  B= Bw+ (-9.9348e-7+ 2.0816e-8*T+ 9.1697e-10*T**2)*S
  bulk_modulus= k0+ A*P+ B*P**2
!
  A= 8.24493e-1- 4.0899e-3*T+ 7.6438e-5*T**2- &
   &   8.2467e-7*T**3+5.3875e-9*T**4
  B= -5.72466e-3 + 1.0227e-4*T- 1.6546e-6*T**2
  C= 4.8314e-4
  rho_w= 999.842594 + 6.793952e-2*T- 9.095290e-3*T**2+ &
   &       1.001685e-4*T**3-1.120083e-6*T**4+ &
   &       6.536336e-9*T**5
  rho_zero= rho_w+ A*S + B*(S**1.5)+ C*(S**2)
!
  RHO= rho_zero/(1- (P/bulk_modulus))
!
END
!
! =========================================================================
!
SUBROUTINE LININT(nx,xtab,ytab,x,y)
!
! Given a value of x return a value of y based on interpolation
! within a table of y values (ytab) corresponding to the x values
! contained in the array xtab.  The subroutine assumes that the
! values in xtab increase monotonically
!
! John Mahaffy 2/12/95
! Modified slightly TRC 2014
!
  integer nx
  double precision xtab(nx), ytab(nx), x, y
!
! local variables
!
  integer i, i1
  double precision  wx
!
  if (x.lt.(xtab(1)).or.x.GT.(xtab(nx))) then
    write(6, *) 'x = ', x, '  is out of table range'
    stop
  endif
  do 100 i=2,nx
       if (x.le.xtab(i)) go to 200
  100 continue
  200 i1=i-1
!
  wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
  y=(1-wx)*ytab(i1)+wx*ytab(i1+1)
!
END
!
! =========================================================================
!
! Dummy routine for ODEPACK. Necessary for Jacobian matrix if stiff ODEs.
!
SUBROUTINE jex()
  RETURN
END
