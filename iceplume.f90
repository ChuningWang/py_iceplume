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
  real(r8) :: Qbar = 250.0d0, Qt = 0.0d0, Qs = 0.0d0
  real(r8) :: dx, dy
!
! Local variable declarations.
!
  real(r8) :: pr, prRef
!
! ==================================================================
! Write header information
! ==================================================================
!
  write(*, *)  ''
  write(*, *)  '==================================================='
  write(*, *)  'ICEPLUME Stand-alone Version'
  write(*, *)  'v1.0'
  write(*, *)  'Chuning Wang'
  write(*, *)  '==================================================='
  write(*, *)  ''
!
! ==================================================================
! Allocate variables
! ==================================================================
!
  CALL allocate_iceplume(ng)
!
! ==================================================================
! Read in some scalar parameters
! ==================================================================
!
  Nr = N(ng)
  NTr = NT(ng)
  dy = 200.0d0
  dx = 200.0d0
!
! ==================================================================
! Read in profiles from OCEAN.
! ==================================================================
!
  write(*, *)  'Read input profiles'
!
! zw
!
  open(unit=5, file='./data/iceplume_zw.txt')
  DO K = 0, Nr
    read(5, *)  PLUME(ng) % zW(I, K)
  ENDDO
  close(unit=5)
!
! temp and salt
!
  open(unit=5, file='./data/iceplume_s.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % sAm(I, K)
  ENDDO
  close(unit=5)
  open(unit=5, file='./data/iceplume_t.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % tpAm(I, K)
  ENDDO
  close(unit=5)
!
! u/v, w
!
  open(unit=5, file='./data/iceplume_v.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % vAm(I, K)
  ENDDO
  close(unit=5)
  open(unit=5, file='./data/iceplume_w.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % wAm(I, K)
  ENDDO
  close(unit=5)
!
! tracers
!
  open(unit=5, file='./data/iceplume_dye01.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % trcAm(I, K, 3)
  ENDDO
  close(unit=5)
!
  DO K = 1, Nr
!
! convert potential temp to in-situ temp
!
    prRef = 101.d3*1.d-4
    pr = prRef + &
      & (abs(PLUME(ng) % zW(I, K))*rho_ref*g)*1.d-4  ! [dbar]
    CALL SW_TEMP(PLUME(ng) % sAm(I, K),  &
               & PLUME(ng) % tpAm(I, K), &
               & pr, prRef, PLUME(ng) % tAm(I, K))
    !PLUME(ng) % tAm(I, K) = PLUME(ng) % tpAm(I, K)
  ENDDO
!
! Discharge tracer concentration
!
  IF (useTracers .and. useInputTracers) THEN
    DO iTracer = 3, NTr
      PLUME(ng) % trcIni(I, iTracer) = 0.0d0
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
  write(*, *)  'Call ICEPLUME main function'
!
  CALL iceplume_calc(ng, I, dx, dy, Qbar, Qt, Qs)
!
! ==================================================================
! Write to file
! ==================================================================
!
  DO K = 1, Nr
    PLUME(ng) % rhoAm(I, K) = RHO(PLUME(ng) % tAm(I, K), &
                                & PLUME(ng) % sAm(I, K), &
                                & PLUME(ng) % zR(I, K)) - 1000.d0
  ENDDO
  DO K = 0, Nr
    PLUME(ng) % rho(I, K) = RHO(PLUME(ng) % t(I, K), &
                              & PLUME(ng) % s(I, K), &
                              & PLUME(ng) % zW(I, K)) - 1000.d0
  ENDDO
  write(*, *)  'Write output to files'
!
  open(unit=15, file='./outputs/plume_out_zw.txt')
  write(15, '(A3, 99 A12)')  'lev', 'zW', &
    & 't', 's', 'r', 'w', 'a', 'mInt', 'volFlux', 'rho'
  DO K = 0, Nr
    write(15, '(I4, 99 E16.8)')  K, PLUME(ng) % zW(I, K), &
      & PLUME(ng) % t(I, K), PLUME(ng) % s(I, K), &
      & PLUME(ng) % r(I, K), PLUME(ng) % w(I, K), &
      & PLUME(ng) % a(I, K), PLUME(ng) % mInt(I, K), &
      & PLUME(ng) % volFlux(I, K), PLUME(ng) % rho(I, K)
  ENDDO
  close(unit=15)
  open(unit=15, file='./outputs/plume_out_zr.txt')
  write(15, '(A3, 99 A12)')  'lev', 'zR', &
    & 'tAm', 'sAm', 'vAm', 'wAm', 'ent', 'det', &
    & 'fwFlux', 'heatFlux', 'rhoAm'
  DO K = 1, Nr
    write(15, '(I4, 99 E16.8)')  K, PLUME(ng) % zR(I, K), &
      & PLUME(ng) % tAm(I, K), PLUME(ng) % sAm(I, K), &
      & PLUME(ng) % vAm(I, K), PLUME(ng) % wAm(I, K), &
      & PLUME(ng) % ent(I, K), PLUME(ng) % det(I, K), &
      & PLUME(ng) % fwFlux(I, K), PLUME(ng) % heatFlux(I, K), &
      & PLUME(ng) % rhoAm(I, K)
  ENDDO
  close(unit=15)
END PROGRAM
