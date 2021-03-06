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
!
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
! tGade     - equavilent temperature of meltwater (Gade Line) [degC]   !
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
! detSigma  - detrainment Gaussian profile standard deviation          !
! avgExp    - logical switch to activate temporal exponential average  !
! avgAlpha  - average weight for exponential smoothing (0 to 1)        !
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
        real(r8), parameter :: tGade      = -(L - cI*tIce)/cW
        real(r8), parameter :: lambda1    = -0.0573_r8
        real(r8), parameter :: lambda2    = 0.0832_r8
        real(r8), parameter :: lambda3    = 0.000761_r8
!
        real(r8), parameter :: GamT       = 0.0220_r8
        real(r8), parameter :: GamS       = 0.000620_r8
        real(r8), parameter :: Cd         = 0.065_r8
! Original value from Cowton et al. 2015
!        real(r8), parameter :: Cd         = 0.0025_r8
!
        real(r8), parameter :: RiB        = 1.0_r8
        real(r8), parameter :: gRedBkg    = 0.01_r8
        real(r8), parameter :: CdBkg      = 0.065_r8
! Original value from Cowton et al. 2015
!        real(r8), parameter :: CdBkg      = 0.0025_r8
        real(r8), parameter :: velBkg     = 0.03_r8
        real(r8), parameter :: wIni       = 0.1_r8
!
        real(r8), parameter :: detSigma   = 0.5_r8
        logical :: avgExp
        real(r8), parameter :: avgAlpha   = 1.0_r8
!
! =====================================================================!
!                                                                      !
! PLUME Type variables                                                 !
!                                                                      !
! Isrc, Jsrc - grid indices acquired from river file.                  !
! dx, dy     - grid lengths in cross/along glacier direction.          !
! dir        - direction of plume. +1 for positve direction and        !
!              -1 for negative direction. 0 for other situation.       !
! trs        - depth integrated total transport [m^3 s^-1]             !
! avgIstr    - Location indices (Istr) when averaging density profiles.!
! avgIend    - Location indices (Iend) when averaging density profiles.!
! avgJstr    - Location indices (Jstr) when averaging density profiles.!
! avgJend    - Location indices (Jend) when averaging density profiles.!
! avgArea    - surface area of average span for each plume source.     !
!                                                                      !
! Profiles                                                             !
!                                                                      !
! zW         - depth [m]                                               !
! f          - plume vertical volume flux [m^3 s^-1]                   !
! w          - plume vertical velocity [m s^-1]                        !
! t          - plume temperature [degC]                                !
! s          - plume salinity [PSU]                                    !
! a          - plume area integrated [m^2]                             !
! mInt       - plume area integrated melt [m^3 s^-1]                   !
! rho        - plume density [kg m^-3]                                 !
!                                                                      !
! zR         - depth at Rho points [m]                                 !
! sAm        - ambient salinity [PSU]                                  !
! tAm        - ambient temperature [degC]                              !
! vAm        - horizontal velocity parallel to glacier [m s^-1]        !
! wAm        - vertical velocity parallel to glacier [m s^-1]          !
! tpAm       - ambient potential temperature [degC]                    !
! rhoAm      - ambient density [kg m^-3]                               !
! rhoAm0     - ambient density of the previous step [kg m^-3]          !
!                                                                      !
! lm         - plume/glacier contact length [m]                        !
! lc         - plume/water contact length [m]                          !
!                                                                      !
! ent        - entrainment rate [m^3 s^-1]                             !
! det        - detrainment rate [m^3 s^-1]                             !
! detI       - detrainment flag                                        !
! detFrac    - fraction of detrainment in vertical direction           !
!                                                                      !
! m          - plume melt rate [m^3 s^-1]                              !
! mB         - background melt rate [m^3 s^-1]                         !
!                                                                      !
! dz         - RHO layer thickness [m]                                 !
!                                                                      !
! Passive tracers                                                      !
!                                                                      !
! trc        - tracer concentration                                    !
! trcAm      - ambient tracer concentration                            !
! trcB       - background meltwater tracer concentration               !
! trcCum     - accumulative tracer concentration                       !
! trcIni     - initial tracer concentration in discharge               !
!                                                                      !
! =====================================================================!
!
        TYPE T_PLUME
!
! Variables.
!
!
! Depth independent variables.
!
          integer,  pointer :: Isrc(:)
          integer,  pointer :: Jsrc(:)
          real(r8), pointer :: dx(:)
          real(r8), pointer :: dy(:)
          real(r8), pointer :: dir(:)
          real(r8), pointer :: trs(:)
          integer,  pointer :: avgIstr(:)
          integer,  pointer :: avgIend(:)
          integer,  pointer :: avgJstr(:)
          integer,  pointer :: avgJend(:)
          real(r8), pointer :: avgArea(:, :)
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
          real(r8), pointer :: rhoAm0(:, :)
!
          real(r8), pointer :: tAm2(:, :)
          real(r8), pointer :: sAm2(:, :)
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
          integer,  pointer :: detI(:, :)
          real(r8), pointer :: detFrac(:, :)
!
! Melt rate, freshwater and heat fluxes.
!
          real(r8), pointer :: m(:, :)
          real(r8), pointer :: mB(:, :)
!
! Other profiles.
!
          real(r8), pointer :: dz(:, :)
!
! Passive tracer concentration.
!
          real(r8), pointer :: trcAm(:, :, :)
          real(r8), pointer :: trcB(:, :)
          real(r8), pointer :: trc(:, :)
          real(r8), pointer :: trcCum(:, :)
          real(r8), pointer :: trcIni(:, :)
        END TYPE T_PLUME
!
        TYPE (T_PLUME), allocatable :: PLUME(:)
!
      CONTAINS
!
      SUBROUTINE allocate_iceplume(ng)
!
! =====================================================================!
!                                                                      !
! Allocate PLUME Type variables.                                       !
!                                                                      !
! =====================================================================!
!
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
!
      integer :: is, itrc, k
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng .EQ. 1) allocate( PLUME(Ngrids) )
!
!  Allocate iceplume Sources/Sinks variables.
!
      allocate ( PLUME(ng) % Isrc (Nsrc(ng)) )
      allocate ( PLUME(ng) % Jsrc (Nsrc(ng)) )
      allocate ( PLUME(ng) % dx   (Nsrc(ng)) )
      allocate ( PLUME(ng) % dy   (Nsrc(ng)) )
      allocate ( PLUME(ng) % dir  (Nsrc(ng)) )
      allocate ( PLUME(ng) % trs  (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgIstr (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgIend (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgJstr (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgJend (Nsrc(ng)) )
      allocate ( PLUME(ng) % avgArea (Nsrc(ng), N(ng)) )
!
      allocate ( PLUME(ng) % zW   (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % f    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % w    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % t    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % s    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % a    (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % mInt (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % rho  (Nsrc(ng), 0:N(ng)) )
!
      allocate ( PLUME(ng) % zR    (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % tAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % sAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % vAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % wAm   (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % tpAm  (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % rhoAm (Nsrc(ng), N(ng)+1) )
      allocate ( PLUME(ng) % rhoAm0 (Nsrc(ng), N(ng)+1) )
!
      allocate ( PLUME(ng) % tAm2  (Nsrc(ng), N(ng)  ) )
      allocate ( PLUME(ng) % sAm2  (Nsrc(ng), N(ng)  ) )
!
      allocate ( PLUME(ng) % lm (Nsrc(ng), 0:N(ng)) )
      allocate ( PLUME(ng) % lc (Nsrc(ng), 0:N(ng)) )
!
      allocate ( PLUME(ng) % ent     (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % det     (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detI    (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % detFrac (Nsrc(ng), N(ng)) )
!
      allocate ( PLUME(ng) % m  (Nsrc(ng), N(ng)) )
      allocate ( PLUME(ng) % mB (Nsrc(ng), N(ng)) )
!
      allocate ( PLUME(ng) % dz (Nsrc(ng), N(ng)) )
!
      allocate ( PLUME(ng) % trcAm  (Nsrc(ng), N(ng), NT(ng)) )
!
      allocate ( PLUME(ng) % trcB   (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trc    (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trcCum (Nsrc(ng), NT(ng)) )
      allocate ( PLUME(ng) % trcIni (Nsrc(ng), NT(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrc(ng)
        PLUME(ng) % Isrc (is) = 0
        PLUME(ng) % Jsrc (is) = 0
        PLUME(ng) % dx   (is) = IniVal
        PLUME(ng) % dy   (is) = IniVal
        PLUME(ng) % dir  (is) = IniVal
        PLUME(ng) % trs  (is) = IniVal
        PLUME(ng) % avgIstr (is) = 0
        PLUME(ng) % avgIend (is) = 0
        PLUME(ng) % avgJstr (is) = 0
        PLUME(ng) % avgJend (is) = 0
        DO k=1,N(ng)
          PLUME(ng) % avgArea (is, k) = IniVal
        ENDDO
!
        DO k=0,N(ng)
          PLUME(ng) % zW   (is, k) = IniVal
          PLUME(ng) % f    (is, k) = IniVal
          PLUME(ng) % w    (is, k) = IniVal
          PLUME(ng) % t    (is, k) = IniVal
          PLUME(ng) % s    (is, k) = IniVal
          PLUME(ng) % a    (is, k) = IniVal
          PLUME(ng) % mInt (is, k) = IniVal
          PLUME(ng) % rho  (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % zR    (is, k) = IniVal
          PLUME(ng) % tAm   (is, k) = IniVal
          PLUME(ng) % sAm   (is, k) = IniVal
          PLUME(ng) % vAm   (is, k) = IniVal
          PLUME(ng) % wAm   (is, k) = IniVal
          PLUME(ng) % tpAm  (is, k) = IniVal
          PLUME(ng) % rhoAm (is, k) = IniVal
          PLUME(ng) % tAm2  (is, k) = IniVal
          PLUME(ng) % sAm2  (is, k) = IniVal
        ENDDO
        PLUME(ng) % rhoAm (is, N(ng)+1) = IniVal
        DO k=1,N(ng)
          PLUME(ng) % rhoAm0 (is, k) = IniVal
        ENDDO
        PLUME(ng) % rhoAm0 (is, N(ng)+1) = IniVal
!
        DO k=0,N(ng)
          PLUME(ng) % lm (is, k) = IniVal
          PLUME(ng) % lc (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % ent     (is, k) = IniVal
          PLUME(ng) % det     (is, k) = IniVal
          PLUME(ng) % detI    (is, k) = 0
          PLUME(ng) % detFrac (is, k) = IniVal
        ENDDO
!
        DO k=1,N(ng)
          PLUME(ng) % m  (is, k) = IniVal
          PLUME(ng) % mB (is, k) = IniVal
!
          PLUME(ng) % dz (is, k) = IniVal
!
          DO itrc=1,NT(ng)
            PLUME(ng) % trcAm (is, k, itrc) = IniVal
          ENDDO
        ENDDO
!
        DO itrc=1,NT(ng)
          PLUME(ng) % trcB   (is, itrc) = IniVal
          PLUME(ng) % trc    (is, itrc) = IniVal
          PLUME(ng) % trcCum (is, itrc) = IniVal
          PLUME(ng) % trcIni (is, itrc) = IniVal
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE allocate_iceplume
      END MODULE mod_iceplume
