PROGRAM iceplume
!
! Load private module
!
  USE mod_iceplume
!
!  Imported variable declarations.
!
  integer :: ng = 1
  integer :: I = 1
  real(r8) :: dx, dy
  real(r8) :: fIni, sIni, tIni, trcIni
  real(r8) :: sgTyp, sgDep, sgLen
  real(r8) :: RHO
!
! Local variable declarations.
!
  real(r8) :: pr, prRef, Nreal
!
! ==================================================================
! Write header information
! ==================================================================
!
  write(*, *)  ''
  write(*, *)  '==================================================='
  write(*, *)  'ICEPLUME Stand-alone Model'
  write(*, *)  'v1.0'
  write(*, *)  'Chuning Wang, chuning@marine.rutgers.edu'
  write(*, *)  '==================================================='
  write(*, *)  ''
!
! ==================================================================
! Read in some scalar parameters
! ==================================================================
!
  open(unit=5, file='./inputs/iceplume_scalar.txt')
  read(5, *)  Nreal, dx, dy, dt(ng),    &
            & fIni, sIni, tIni, trcIni, &
            & sgTyp, sgDep, sgLen
  close(unit=5)
  N(ng) = INT(nreal)
  Nr = N(ng)
  NTr = NT(ng)
!
! ==================================================================
! Allocate variables
! ==================================================================
!
  CALL allocate_iceplume(ng)
!
! ==================================================================
! Read in profiles from OCEAN.
! ==================================================================
!
  write(*, *)  'Reading input data...'
!
! zw
!
  open(unit=5, file='./inputs/iceplume_zw.txt')
  DO K = 0, Nr
    read(5, *)  PLUME(ng) % zW(I, K)
  ENDDO
  close(unit=5)
!
! zr
!
  open(unit=5, file='./inputs/iceplume_zr.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % tpAm(I, K), PLUME(ng) % sAm(I, K),  &
              & PLUME(ng) % vAm(I, K), PLUME(ng) % wAm(I, K),   &
              & PLUME(ng) % trcAm(I, K, 3)
  ENDDO
  close(unit=5)
!
! convert potential temp to in-situ temp
!
  DO K = 1, Nr
    IF (usePotTemp) THEN
      prRef = 101.d3*1.d-4
      pr = prRef + &
        & (abs(PLUME(ng) % zW(I, K))*rhoRef*g)*1.d-4  ! [dbar]
      CALL SW_TEMP(PLUME(ng) % sAm(I, K),  &
                 & PLUME(ng) % tpAm(I, K), &
                 & pr, prRef, PLUME(ng) % tAm(I, K))
    ELSE
      PLUME(ng) % tAm(I, K) = PLUME(ng) % tpAm(I, K)
    ENDIF
  ENDDO
!
! Discharge tracer concentration
!
  IF (useTracers .and. useInputTracers) THEN
    DO iTracer = 3, NTr
      PLUME(ng) % trcIni(I, iTracer) = trcIni
    ENDDO
  ELSE
    DO iTracer = 3, NTr
      PLUME(ng) % trcIni(I, iTracer) = 0.0d0
    ENDDO
  ENDIF
!
! ==================================================================
! Call main function
! ==================================================================
!
  write(*, *)  'Calculating ICEPLUME...'
!
  CALL ICEPLUME_CALC(ng, I, dx, dy,       &
                   & fIni, tIni, sIni,    &
                   & INT(sgTyp), sgDep, sgLen)
!
! ==================================================================
! Write to file
! ==================================================================
!
  DO K = 1, Nr
    PLUME(ng) % rhoAm(I, K) = RHO(PLUME(ng) % tAm(I, K), &
                                & PLUME(ng) % sAm(I, K), &
                                & PLUME(ng) % zR(I, K)) - 1000.0d0
  ENDDO
  DO K = 0, Nr
    PLUME(ng) % rho(I, K) = RHO(PLUME(ng) % t(I, K), &
                              & PLUME(ng) % s(I, K), &
                              & PLUME(ng) % zW(I, K)) - 1000.0d0
  ENDDO
  write(*, *)  'Writing output to files...'
!
  open(unit=15, file='./outputs/iceplume_zw.txt')
  write(15, '(A4, 99 A20)')  'lev', 'zW',                 &
    & 'f', 'w', 't', 's', 'a', 'mInt', 'rho'
  DO K = 0, Nr
    write(15, '(I4, 99 E20.8)')  K, PLUME(ng) % zW(I, K), &
      & PLUME(ng) % f(I, K), PLUME(ng) % w(I, K),         &
      & PLUME(ng) % t(I, K), PLUME(ng) % s(I, K),         &
      & PLUME(ng) % a(I, K), PLUME(ng) % mInt(I, K),      &
      & PLUME(ng) % rho(I, K)
  ENDDO
  close(unit=15)
  open(unit=15, file='./outputs/iceplume_zr.txt')
  write(15, '(A4, 99 A20)')  'lev', 'zR', 'dz',           &
    & 'ent', 'det', 'tAm', 'sAm', 'vAm', 'wAm',           &
    & 'rhoAm', 'm', 'mB'
  DO K = 1, Nr
    write(15, '(I4, 99 E20.8)')                           &
      & K, PLUME(ng) % zR(I, K), PLUME(ng) % dz(I, K),    &
      & PLUME(ng) % ent(I, K), PLUME(ng) % det(I, K),     &
      & PLUME(ng) % tAm(I, K), PLUME(ng) % sAm(I, K),     &
      & PLUME(ng) % vAm(I, K), PLUME(ng) % wAm(I, K),     &
      & PLUME(ng) % rhoAm(I, K), PLUME(ng) % m(I, K),     &
      & PLUME(ng) % mB(I, K)
  ENDDO
  close(unit=15)
  open(unit=15, file='./outputs/iceplume_dye.txt')
    write(15, '(99 E20.8)') PLUME(ng) % trc(I, :)
  close(unit=15)
END PROGRAM
