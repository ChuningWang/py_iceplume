MODULE mod_iceplume
! This module stores all global variables.
  USE mod_kinds
  USE mod_iceplume_py
  implicit none

! ====================================================================
! User defined variables and logical switches
! These parameters are read when coupling with GCMs.
  logical :: usePlumeDiagnostics = .false.
  logical :: useConePlume        = .true.
  logical :: useSheetPlume       = .false.
  logical :: conserveMass        = .false.
  logical :: correctMeltEnt      = .false.
  logical :: useBkgMelt          = .false.
  logical :: useTracers          = .true.
  logical :: useInputTracers     = .true.
  logical :: useDebug            = .true.
  logical :: checkCFL            = .true.
  logical :: checkCFLVerbal      = .false.
  logical :: checkRiB            = .true.
  integer :: depthFinder         = 2

  integer :: ngr              ! nested grid ID
  integer :: Nr               ! total number of layers
  integer :: NTr              ! total number of tracers
  real(r8) :: dy              ! grid cell resolution [m]
  real(r8) :: dx              ! grid cell resolution [m]
  real(r8) :: dtr             ! slow-mode time step [m],
                              ! this is required when
                              ! checkCFL is activated.
  real(r8) :: CuMax = 0.8d0   ! maximum Cu number to triger ckeckCFL.
  real(r8) :: RiBMin = 0.5d0  ! minimum Richardson number to triger
                              ! checkRiB.
! ice bottom depth [m], if iceDepth is positive,
! then the ice bottom depth is always equal to
! the water depth (iceDepthK = 1)
!  real(r8) :: iceDepth = -200.d0
  real(r8) :: iceDepth = 1.d0

! ====================================================================

! Model parameters
  real(r8), parameter :: pi            = 4.0d0*atan(1.0d0)    ! Pi

! For plume model
! -------------- Parameters ------------------------------------------
! E_0 - entrainment rate (alpha)
! iceTemp - ice temperature [degC]
! rho_ref - reference density [kg m^-3]
! rho_fresh - reference density [kg m^-3]
! rho_ice - ice density [kg m^-3]
! g - gravity acceleration [m s^-2]
! c_w - heat capacity of water [J kg^-1 degC^-1]
! c_i - heat capacity of ice [J kg^-1 degC^-1]
! L - latent heat of melting [J kg^-1]
! lambda1 - freezing point slope [degC PSU^-1]
! lambda2 - freezing point offset [degC]
! lambda3 - freezing point depth slope [degC m^-1]
! GamT - thermal turbulent transfer coefficient
! GamS - salt turbulent transfer coefficient
! Cd - ice-plume drag coefficient
! BackgroundVel - background velocity [m s^-1]
!
! wIni - Initial (discharge) velocity [m s^-1]
!
! Meltwater entrainment ratio [ratio]
! This is a parameter to quantify how much background meltwater is
! entrained into the plume. This ratio varies between 0 and 1, where
! 0 means no entrainment into the plume and 1 means all entrainment
!
! iceDepthK - ice bottom layer index [sigma points]
! plumeDepthK - neutral buoyancy plume layer index [sigma points]
! mass2rUnit - reciprocal density [kg^-1 m^3]
! -------------- Parameters ------------------------------------------

  real(r8), parameter :: E_0           = 1.d-1
  real(r8), parameter :: iceTemp       = -1.d1
  real(r8), parameter :: rho_ref       = 1.020d3
  real(r8), parameter :: rho_fresh     = 1.000d3
  real(r8), parameter :: rho_ice       = 0.9167d3
  real(r8), parameter :: g             = 9.81d0
  real(r8), parameter :: c_w           = 3.974d3
  real(r8), parameter :: c_i           = 2.009d3
  real(r8), parameter :: L             = 3.35d5
  real(r8), parameter :: lambda1       = -5.73d-2
  real(r8), parameter :: lambda2       = 8.32d-2
  real(r8), parameter :: lambda3       = 7.61d-4

  real(r8), parameter :: GamT          = 2.20d-2
  real(r8), parameter :: GamS          = 6.20d-4
  real(r8), parameter :: Cd            = 2.50d-3
  real(r8), parameter :: backgroundVel = 1.d-2

  real(r8), parameter :: wIni          = 1.d0

  real(r8), parameter :: meltEnt       = 5.d-1

  integer :: iceDepthK                             
  integer :: plumeDepthK                           

! Taken from pkg/icefront
! real(r8), parameter :: mass2rUnit = 1.d0/rho_ref

! Looping vars
  integer :: K  ! vertical layers
  integer :: iTracer  ! tracer index

! ====================================================================
! Profiles
! ====================================================================
! AT RHO POINTS
! sAm           - ambient salinity [PSU]
! tAm           - ambient temperature [degC]
! vAm           - horizontal velocity parallel to the glacier wall [m s^-1]
! wAm           - vertical velocity parallel to the glacier wall [m s^-1]
! tpAm          - ambient potential temperature [degC]
! zRho          - depth at Rho points [m]
! pr            - pressure [dbar]
! mAv           - total melt rate [m s^-1]
! m             - plume melt rate [m s^-1]
! mAm           - background melt rate [m s^-1]
! fwFlux        - meltwater freshwater flux into grid cell [kg m^-2 s^-1]
!
! heatFlux      - meltwater heat flux into grid cell [W m^-2]
! tendT         - a tendency term calculated as in MITgcm.
! tendS         - a tendency term calculated as in MITgcm.
! dz            - veritcal resolution [m]
! ent           - entrainment rate [m^3 s^-1]
! det           - detrainment rate [m^3 s^-1]
! detI          - detrainment flag
! bvf           - Brunt-Vaisala frequency [s^-2]
!
! AT OMEGA POINTS
! zW            - depth [m]
! r             - radius of plume [m]
! w             - vertical velocity of plume [m s^-1]
! t             - temperature of plume [degC]
! s             - salinity of plume [PSU]
! a             - area of plume [m^2]
! mInt          - area integrated melt of plume [m^3 s^-1]
! volFLux       - vertical volume flux [m^3 s^-1]
!
! PASSIVE TRACERS
! ptrProf       - passive tracer concentration
! ====================================================================
  TYPE T_PLUME
! plume model state.
! Grid, ambient state.
    real(r8), pointer :: zW(:)
    real(r8), pointer :: tAm(:)
    real(r8), pointer :: sAm(:)
    real(r8), pointer :: vAm(:)
    real(r8), pointer :: wAm(:)
    real(r8), pointer :: tpAm(:)
    
! Plume state.
    real(r8), pointer :: r(:)
    real(r8), pointer :: t(:)
    real(r8), pointer :: s(:)
    real(r8), pointer :: w(:)
    real(r8), pointer :: a(:)
    real(r8), pointer :: mInt(:)

! Melt rate and fluxes.
    real(r8), pointer :: mAv(:)
    real(r8), pointer :: m(:)
    real(r8), pointer :: mAm(:)
    real(r8), pointer :: ent(:)
    real(r8), pointer :: det(:)
    integer(r8), pointer :: detI(:)
    real(r8), pointer :: fwFlux(:)
    real(r8), pointer :: heatFlux(:)

! For Richardson number calculation.
!    real(r8), pointer :: bvf(:)

! Tendency terms.
    real(r8), pointer :: tendT(:)
    real(r8), pointer :: tendS(:)

! Other profiles.
    real(r8), pointer :: zRho(:)
    real(r8), pointer :: dz(:)
    real(r8), pointer :: volFlux(:)
    real(r8), pointer :: tB(:)
    real(r8), pointer :: sB(:)

! Passive tracer concentration.
    real(r8), pointer :: trcAm(:, :)
    real(r8), pointer :: trc(:)
    real(r8), pointer :: trcCum(:)
    real(r8), pointer :: trcIni(:)
  END TYPE T_PLUME

  TYPE (T_PLUME), allocatable :: PLUME(:)

! ===================================================================
  CONTAINS
    SUBROUTINE allocate_iceplume(ng)
      integer :: ng
      IF (ng .EQ. 1) allocate(PLUME(Ngrids))
! Allocate profiles
        allocate(PLUME(ng) % zW(N(ng)+1))
        allocate(PLUME(ng) % tAm(N(ng)))
        allocate(PLUME(ng) % sAm(N(ng)))
        allocate(PLUME(ng) % vAm(N(ng)))
        allocate(PLUME(ng) % wAm(N(ng)))
        allocate(PLUME(ng) % tpAm(N(ng)))
!
        allocate(PLUME(ng) % r(N(ng)+1))
        allocate(PLUME(ng) % t(N(ng)+1))
        allocate(PLUME(ng) % s(N(ng)+1))
        allocate(PLUME(ng) % w(N(ng)+1))
        allocate(PLUME(ng) % a(N(ng)+1))
        allocate(PLUME(ng) % mInt(N(ng)+1))
!
        allocate(PLUME(ng) % mAv(N(ng)))
        allocate(PLUME(ng) % m(N(ng)))
        allocate(PLUME(ng) % mAm(N(ng)))
        allocate(PLUME(ng) % ent(N(ng)))
        allocate(PLUME(ng) % det(N(ng)))
        allocate(PLUME(ng) % detI(N(ng)))
        allocate(PLUME(ng) % fwFlux(N(ng)))
        allocate(PLUME(ng) % heatFlux(N(ng)))
!
!        allocate(PLUME(ng) % bvf(N(ng)))
!
        allocate(PLUME(ng) % tendT(N(ng)))
        allocate(PLUME(ng) % tendS(N(ng)))
!
        allocate(PLUME(ng) % zRho(N(ng)))
        allocate(PLUME(ng) % dz(N(ng)))
        allocate(PLUME(ng) % volFlux(N(ng)+1))
        allocate(PLUME(ng) % tB(N(ng)))
        allocate(PLUME(ng) % sB(N(ng)))
!
        allocate(PLUME(ng) % trcAm(N(ng), NT(ng)))
        allocate(PLUME(ng) % trc(NT(ng)))
        allocate(PLUME(ng) % trcCum(NT(ng)))
        allocate(PLUME(ng) % trcIni(NT(ng)))
    END SUBROUTINE allocate_iceplume
END
