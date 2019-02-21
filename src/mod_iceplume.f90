MODULE mod_iceplume
! ==================================================================!
!                                                                   !
! These are the module functions of  model.                 !
!                                                                   !
! ==================================================================!
!
! This module stores all global variables.
!
  USE mod_kinds
  USE mod_py_iceplume
  implicit none
!
! ====================================================================!
!                                                                     !
! User defined variables and logical switches                         !
! These parameters are read when coupling with GCMs.                  !
!                                                                     !
! ====================================================================!
!
  logical :: useTracers          = .true.
  logical :: useInputTracers     = .true.
  logical :: useBkgMelt          = .true.
  logical :: useMeltTracers      = .true.
  logical :: checkCFL            = .true.
  logical :: checkRiB            = .true.
  integer :: depthFinder         = 2
!
! ====================================================================!
!                                                                     !
! Model parameters                                                    !
!                                                                     !
! alpha     - entrainment rate                                        !
! tIce      - ice temperature [degC]                                  !
! sIce      - ice salinity [PSU]                                      !
! rhoRef    - reference density [kg m^-3]                             !
! rhoAir    - air density [kg m^-3]                                   !
! g         - gravity acceleration [m s^-2]                           !
! cW        - heat capacity of water [J kg^-1 degC^-1]                !
! cI        - heat capacity of ice [J kg^-1 degC^-1]                  !
! L         - latent heat of melting [J kg^-1]                        !
! lambda1   - freezing point slope [degC PSU^-1]                      !
! lambda2   - freezing point offset [degC]                            !
! lambda3   - freezing point depth slope [degC m^-1]                  !
! GamT      - thermal turbulent transfer coefficient                  !
! GamS      - salt turbulent transfer coefficient                     !
! Cd        - ice-plume drag coefficient                              !
!                                                                     !
! N2Bkg     - background N2                                           !
! CdBkg     - background ice-plume drag coefficient                   !
! velBkg    - background velocity [m s^-1]                            !
!                                                                     !
! wIni      - initial (discharge) velocity [m s^-1]                   !
!                                                                     !
! CuMax     - maximum Courant number to triger checkCFL.              !
! RiBMin    - minimum Richardson number to triger checkRiB.           !
!                                                                     !
! iceDepth  - ice bottom depth [m]                                    !
!             if iceDepth is positive (iceDepth = 1),                 !
!             then iceDepth = water depth                             !
! L0        - finite line source length [m]                           !
!                                                                     !
! ====================================================================!
!
  real(r8), parameter :: pi = 4.0d0*ATAN(1.0d0)    ! Pi
!
  real(r8), parameter :: alpha      = 1.d-1
  real(r8), parameter :: tIce       = -1.d1
  real(r8), parameter :: sIce       = 0.d0
  real(r8), parameter :: rhoRef     = 1.020d3
  real(r8), parameter :: rhoAir     = 1.225
  real(r8), parameter :: g          = 9.81d0
  real(r8), parameter :: cW         = 3.974d3
  real(r8), parameter :: cI         = 2.d3
  real(r8), parameter :: L          = 3.35d5
  real(r8), parameter :: lambda1    = -5.73d-2
  real(r8), parameter :: lambda2    = 8.32d-2
  real(r8), parameter :: lambda3    = 7.61d-4
!
  real(r8), parameter :: GamT       = 2.20d-2
  real(r8), parameter :: GamS       = 6.20d-4
  real(r8), parameter :: Cd         = 6.5d-2
!
  real(r8), parameter :: N2Bkg      = 1.d-3
  real(r8), parameter :: CdBkg      = 2.5d-3
  real(r8), parameter :: velBkg     = 1.d-2
!
  real(r8), parameter :: wIni       = 1.d0
!
  real(r8), parameter :: CuMax      = 0.5d0
  real(r8), parameter :: RiBMin     = 1.0d0
!
! ====================================================================!
!                                                                     !
! Variables                                                           !
!                                                                     !
! dir           - direction of plume. +1 for positve direction and    !
!                 -1 for negative direction. 0 for other situation.   !
!                                                                     !
! Profiles                                                            !
!                                                                     !
! zR            - depth at Rho points [m]                             !
! sAm           - ambient salinity [PSU]                              !
! tAm           - ambient temperature [degC]                          !
! vAm           - horizontal velocity parallel to glacier [m s^-1]    !
! wAm           - vertical velocity parallel to glacier [m s^-1]      !
! tpAm          - ambient potential temperature [degC]                !
! rhoAm         - ambient density [kg m^-3]                           !
!                                                                     !
! zW            - depth [m]                                           !
! f             - plume vertical volume flux [m^3 s^-1]               !
! w             - plume vertical velocity [m s^-1]                    !
! t             - plume temperature [degC]                            !
! s             - plume salinity [PSU]                                !
! a             - plume area integrated [m^2]                         !
! mInt          - plume area integrated melt [m^3 s^-1]               !
!                                                                     !
! lm            - plume/glacier contact length [m]                    !
! lc            - plume/water contact length [m]                      !
!                                                                     !
! ent           - entrainment rate [m^3 s^-1]                         !
! det           - detrainment rate [m^3 s^-1]                         !
! detI          - detrainment flag                                    !
! verI          - vertical acceleration flag                          !
!                                                                     !
! m             - plume melt rate [m^3 s^-1]                          !
! mB            - background melt rate [m^3 s^-1]                     !
!                                                                     !
! dz            - RHO layer thickness [m]                             !
!                                                                     !
! Passive tracers                                                     !
! trc           - passive tracer concentration                        !
! trcAm         - ambient passive tracer concentration                !
! trcB          - meltwater passive tracer concentration              !
! trcCum        - accumulative passive tracer concentration           !
! trcIni        - initial passive tracer concentration in discharge   !
!                                                                     !
! trcAmToB      - tracer flux rate associated with background         !
!               - melt [unit s^-1]                                    !
!                                                                     !
! ====================================================================!
!
  TYPE T_PLUME
!
! Variables.
!
    real(r8), pointer :: dir(:)
!
! Grid, ambient state.
!
    real(r8), pointer :: zR(:, :)
    real(r8), pointer :: tAm(:, :)
    real(r8), pointer :: sAm(:, :)
    real(r8), pointer :: vAm(:, :)
    real(r8), pointer :: wAm(:, :)
    real(r8), pointer :: tpAm(:, :)
    real(r8), pointer :: rhoAm(:, :)
!
! Plume state.
!
    real(r8), pointer :: zW(:, :)
    real(r8), pointer :: f(:, :)
    real(r8), pointer :: w(:, :)
    real(r8), pointer :: t(:, :)
    real(r8), pointer :: s(:, :)
    real(r8), pointer :: a(:, :)
    real(r8), pointer :: mInt(:, :)
    real(r8), pointer :: rho(:, :)
!
! Plume shape parameters.
!
    real(r8), pointer :: lm(:, :)
    real(r8), pointer :: lc(:, :)
!
! Volume fluxes.
!
    real(r8), pointer :: ent(:, :)
    real(r8), pointer :: det(:, :)
    integer(r8), pointer :: detI(:, :)
    integer(r8), pointer :: verI(:, :)
!
! Melt rate, freshwater and heat fluxes.
!
    real(r8), pointer :: m(:, :)
    real(r8), pointer :: mB(:, :)
    real(r8), pointer :: trcAmToB(:, :, :)
!
! Other profiles.
!
    real(r8), pointer :: dz(:, :)
!
! Passive tracer concentration.
!
    real(r8), pointer :: trcAm(:, :, :)
    real(r8), pointer :: trcB(:, :, :)
    real(r8), pointer :: trc(:, :)
    real(r8), pointer :: trcCum(:, :)
    real(r8), pointer :: trcIni(:, :)
  END TYPE T_PLUME
!
  TYPE (T_PLUME), allocatable :: PLUME(:)
!
! ===================================================================
!
  CONTAINS
    SUBROUTINE allocate_iceplume(ng)
      integer :: ng
      IF (ng .EQ. 1) allocate( PLUME(Ngrids) )
!
! Allocate profiles
!
        allocate( PLUME(ng) % dir(1:Nsrc(ng)) )
!
        allocate( PLUME(ng) % zR    (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % tAm   (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % sAm   (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % vAm   (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % wAm   (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % tpAm  (Nsrc(ng), N(ng)  ) )
        allocate( PLUME(ng) % rhoAm (Nsrc(ng), N(ng)+1) )
!
        allocate( PLUME(ng) % zW   (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % f    (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % w    (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % t    (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % s    (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % a    (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % mInt (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % rho  (Nsrc(ng), 0:N(ng)) )
!
        allocate( PLUME(ng) % lm (Nsrc(ng), 0:N(ng)) )
        allocate( PLUME(ng) % lc (Nsrc(ng), 0:N(ng)) )
!
        allocate( PLUME(ng) % ent  (Nsrc(ng), N(ng)) )
        allocate( PLUME(ng) % det  (Nsrc(ng), N(ng)) )
        allocate( PLUME(ng) % detI (Nsrc(ng), N(ng)) )
        allocate( PLUME(ng) % verI (Nsrc(ng), N(ng)) )
!
        allocate( PLUME(ng) % m        (Nsrc(ng), N(ng)) )
        allocate( PLUME(ng) % mB       (Nsrc(ng), N(ng)) )
        allocate( PLUME(ng) % trcAmToB (Nsrc(ng), N(ng), NT(ng)) )
!
        allocate( PLUME(ng) % dz (Nsrc(ng), N(ng)) )
!
        allocate( PLUME(ng) % trcAm  (Nsrc(ng), N(ng), NT(ng)) )
        allocate( PLUME(ng) % trcB   (Nsrc(ng), N(ng), NT(ng)) )
        allocate( PLUME(ng) % trc    (Nsrc(ng), NT(ng)) )
        allocate( PLUME(ng) % trcCum (Nsrc(ng), NT(ng)) )
        allocate( PLUME(ng) % trcIni (Nsrc(ng), NT(ng)) )
    END SUBROUTINE allocate_iceplume
END
