MODULE mod_iceplume
! =====================================================================!
!                                                                      !
! These are the module functions of iceplume model.                    !
!                                                                      !
! =====================================================================!
!
! This module stores all global variables.
!
  USE mod_kinds
  USE mod_py_iceplume
  implicit none
!
! =====================================================================!
!                                                                      !
! Model parameters                                                     !
!                                                                      !
! alpha     - entrainment rate                                         !
! tIce      - ice temperature [degC]                                   !
! sIce      - ice salinity [PSU]                                       !
! rhoRef    - reference density [kg m^-3]                              !
! rhoAir    - air density [kg m^-3]                                    !
! g         - gravity acceleration [m s^-2]                            !
! cW        - heat capacity of water [J kg^-1 degC^-1]                 !
! cI        - heat capacity of ice [J kg^-1 degC^-1]                   !
! L         - latent heat of melting [J kg^-1]                         !
! lambda1   - freezing point slope [degC PSU^-1]                       !
! lambda2   - freezing point offset [degC]                             !
! lambda3   - freezing point depth slope [degC m^-1]                   !
!                                                                      !
! GamT      - thermal turbulent transfer coefficient                   !
! GamS      - salt turbulent transfer coefficient                      !
! Cd        - ice-plume drag coefficient                               !
!                                                                      !
! RiB       - critical Richardson number                               !
! gRedBkg   - background reduced gravity                               !
! CdBkg     - background ice-plume drag coefficient                    !
! velBkg    - background velocity [m s^-1]                             !
! wIni      - initial (discharge) velocity [m s^-1]                    !
!                                                                      !
! =====================================================================!
!
  real(r8), parameter :: pi = 4.0d0*ATAN(1.0d0)    ! Pi
!
  real(r8), parameter :: alpha      = 0.1_r8
  real(r8), parameter :: tIce       = -10.0_r8
  real(r8), parameter :: sIce       = 0.0_r8
  real(r8), parameter :: rhoRef     = 1020.0_r8
  real(r8), parameter :: rhoAir     = 1.225_r8
  real(r8), parameter :: g          = 9.81_r8
  real(r8), parameter :: cW         = 3974.0_r8
  real(r8), parameter :: cI         = 2000.0_r8
  real(r8), parameter :: L          = 335000.0_r8
  real(r8), parameter :: lambda1    = -0.0573_r8
  real(r8), parameter :: lambda2    = 0.0832_r8
  real(r8), parameter :: lambda3    = 0.000761_r8
!
  real(r8), parameter :: GamT       = 0.0220_r8
  real(r8), parameter :: GamS       = 0.000620_r8
  real(r8), parameter :: Cd         = 0.065_r8
!
  real(r8), parameter :: RiB        = 1.0_r8
  real(r8), parameter :: gRedBkg    = 0.1_r8
  real(r8), parameter :: CdBkg      = 0.0025_r8
  real(r8), parameter :: velBkg     = 0.03_r8
  real(r8), parameter :: wIni       = 1.0_r8
!
! =====================================================================!
!                                                                      !
! PLUME Type variables                                                 !
!                                                                      !
! dir           - direction of plume. +1 for positve direction and     !
!                 -1 for negative direction. 0 for other situation.    !
!                                                                      !
! Profiles                                                             !
!                                                                      !
! zW            - depth [m]                                            !
! f             - plume vertical volume flux [m^3 s^-1]                !
! w             - plume vertical velocity [m s^-1]                     !
! t             - plume temperature [degC]                             !
! s             - plume salinity [PSU]                                 !
! a             - plume area integrated [m^2]                          !
! mInt          - plume area integrated melt [m^3 s^-1]                !
! rho           - plume density [kg m^-3]                              !
!                                                                      !
! zR            - depth at Rho points [m]                              !
! sAm           - ambient salinity [PSU]                               !
! tAm           - ambient temperature [degC]                           !
! vAm           - horizontal velocity parallel to glacier [m s^-1]     !
! wAm           - vertical velocity parallel to glacier [m s^-1]       !
! tpAm          - ambient potential temperature [degC]                 !
! rhoAm         - ambient density [kg m^-3]                            !
!                                                                      !
! lm            - plume/glacier contact length [m]                     !
! lc            - plume/water contact length [m]                       !
!                                                                      !
! ent           - entrainment rate [m^3 s^-1]                          !
! det           - detrainment rate [m^3 s^-1]                          !
! detI          - detrainment flag                                     !
! detFrac       - fraction of detrainment in vertical direction        !
!                                                                      !
! m             - plume melt rate [m^3 s^-1]                           !
! mB            - background melt rate [m^3 s^-1]                      !
!                                                                      !
! dz            - RHO layer thickness [m]                              !
!                                                                      !
! Passive tracers                                                      !
!                                                                      !
! trc           - tracer concentration                                 !
! trcAm         - ambient tracer concentration                         !
! trcB          - background meltwater tracer concentration            !
! trcCum        - accumulative tracer concentration                    !
! trcIni        - initial tracer concentration in discharge            !
!                                                                      !
! trcAmToB      - tracer flux rate associated with background          !
!               - melt [unit s^-1]                                     !
!                                                                      !
! wAvg          - average weight for the current time step (0 to 1)    !
!                                                                      !
! =====================================================================!
!
  TYPE T_PLUME
!
! Variables.
!
    real(r8), pointer :: dir(:)
!
! Plume state (omega surface).
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
! Grid, ambient state (rho surface).
!
    real(r8), pointer :: zR(:, :)
    real(r8), pointer :: tAm(:, :)
    real(r8), pointer :: sAm(:, :)
    real(r8), pointer :: vAm(:, :)
    real(r8), pointer :: wAm(:, :)
    real(r8), pointer :: tpAm(:, :)
    real(r8), pointer :: rhoAm(:, :)
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
    real(r8), pointer :: detFrac(:, :)
!
! For the neutral buoyancy detrainment model.
!
    real(r8), pointer :: detF(:, :)
    real(r8), pointer :: detE(:, :)
    real(r8), pointer :: detTrc(:, :, :)
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
!
! For temporally average background density profile.
!
    real(r8) :: wAvg
  END TYPE T_PLUME
!
  TYPE (T_PLUME), allocatable :: PLUME(:)
!
! =====================================================================!
!                                                                      !
! Allocate PLUME Type variables.                                       !
!                                                                      !
! =====================================================================!
!
  CONTAINS
    SUBROUTINE allocate_iceplume(ng)
      integer :: ng, K, is, itrc
      IF (ng .EQ. 1) allocate( PLUME(Ngrids) )
!
! Allocate profiles
!
      allocate( PLUME(ng) % dir(1:Nsrc(ng)) )
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
      allocate( PLUME(ng) % zR    (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % tAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % sAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % vAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % wAm   (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % tpAm  (Nsrc(ng), N(ng)  ) )
      allocate( PLUME(ng) % rhoAm (Nsrc(ng), N(ng)+1) )
!
      allocate( PLUME(ng) % lm (Nsrc(ng), 0:N(ng)) )
      allocate( PLUME(ng) % lc (Nsrc(ng), 0:N(ng)) )
!
      allocate( PLUME(ng) % ent  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % det  (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detI (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detFrac (Nsrc(ng), N(ng)) )
!
      allocate( PLUME(ng) % detF    (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detE    (Nsrc(ng), N(ng)) )
      allocate( PLUME(ng) % detTrc  (Nsrc(ng), N(ng), NT(ng)) )
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
!
      PLUME(ng) % wAvg = 1.0
      DO is = 1, Nsrc(ng)
        DO K = 1, N(ng)
          PLUME(ng) % detF(is, K) = 0.0
          PLUME(ng) % detE(is, K) = 0.0
          DO itrc = 1, NT(ng)
            PLUME(ng) % detTrc(is, K, itrc) = 0.0
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE allocate_iceplume
END
