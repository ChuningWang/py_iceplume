PROGRAM iceplume
!
! Load private module
!
  USE mod_iceplume
!
!  Imported variable declarations.
!
  integer :: ng = 1
  real(r8) :: Qbar = 250.0d0, Qt = 0.0d0, Qs = 0.0d0
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
  ngr = ng
  Nr = 40
  NTr = 3
  dy = 200.0d0
  dx = 200.0d0
!
! If checkCFL is activited, need to read in slow-mode time step
!
  IF (checkCFL) THEN
    dtr = 30.0d0
  ENDIF
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
    read(5, *)  PLUME(ng) % zW(K)
  ENDDO
  close(unit=5)
!
! temp and salt
!
  open(unit=5, file='./data/iceplume_s.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % sAm(K)
  ENDDO
  close(unit=5)
  open(unit=5, file='./data/iceplume_t.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % tpAm(K)
  ENDDO
  close(unit=5)
!
! u/v, w
!
  open(unit=5, file='./data/iceplume_v.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % vAm(K)
  ENDDO
  close(unit=5)
  open(unit=5, file='./data/iceplume_w.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % wAm(K)
  ENDDO
  close(unit=5)
!
! tracers
!
  open(unit=5, file='./data/iceplume_dye01.txt')
  DO K = 1, Nr
    read(5, *)  PLUME(ng) % trcAm(K, 3)
  ENDDO
  close(unit=5)
!
  DO K = 1, Nr
!
! convert potential temp to in-situ temp
!
    prRef = 101.d3*1.d-4
    pr = prRef + &
      & (abs(PLUME(ng) % zW(K))*rho_ref*g)*1.d-4  ! [dbar]
    CALL SW_TEMP(PLUME(ng) % sAm(K),  &
               & PLUME(ng) % tpAm(K), &
               & pr, prRef, PLUME(ng) % tAm(K))
  ENDDO
!
! Discharge tracer concentration
!
  IF (useTracers .and. useInputTracers) THEN
    DO iTracer = 3, NTr
      PLUME(ng) % trcIni(iTracer) = 0.0d0
    ENDDO
  ELSE
    DO iTracer = 3, NTr
      PLUME(ng) % trcIni(iTracer) = 0.d0
    ENDDO
  ENDIF
!
! ==================================================================
! Call main function
! ==================================================================
!
  write(*, *)  'Call ICEPLUME main function'
!
  CALL iceplume_calc(ng, Qbar, Qt, Qs)
!
! ==================================================================
! Write to file
! ==================================================================
!
  DO K = 1, Nr
    PLUME(ng) % rhoAm(K) = RHO(PLUME(ng) % tAm(K), &
                             & PLUME(ng) % sAm(K), &
                             & PLUME(ng) % zR(K)) - 1000.d0
  ENDDO
  DO K = 0, Nr
    PLUME(ng) % rho(K) = RHO(PLUME(ng) % t(K), &
                           & PLUME(ng) % s(K), &
                           & PLUME(ng) % zW(K)) - 1000.d0
  ENDDO
  write(*, *)  'Write output to files'
!
  open(unit=15, file='./outputs/plume_out_zw.txt')
  write(15, '(A3, 99 A12)')  'lev', 'zW', &
    & 't', 's', 'r', 'w', 'a', 'mInt', 'volFlux', 'rho'
  DO K = 0, Nr
    write(15, '(I4, 99 E12.4)')  K, PLUME(ng) % zW(K), &
      & PLUME(ng) % t(K), PLUME(ng) % s(K), &
      & PLUME(ng) % r(K), PLUME(ng) % w(K), &
      & PLUME(ng) % a(K), PLUME(ng) % mInt(K), &
      & PLUME(ng) % volFlux(K), PLUME(ng) % rho(K)
  ENDDO
  close(unit=15)
  open(unit=15, file='./outputs/plume_out_zr.txt')
  write(15, '(A3, 99 A12)')  'lev', 'zR', &
    & 'tAm', 'sAm', 'vAm', 'wAm', 'ent', 'det', &
    & 'fwFlux', 'heatFlux', 'rhoAm'
  DO K = 1, Nr
    write(15, '(I4, 99 E12.4)')  K, PLUME(ng) % zR(K), &
      & PLUME(ng) % tAm(K), PLUME(ng) % sAm(K), &
      & PLUME(ng) % vAm(K), PLUME(ng) % wAm(K), &
      & PLUME(ng) % ent(K), PLUME(ng) % det(K), &
      & PLUME(ng) % fwFlux(K), PLUME(ng) % heatFlux(K), &
      & PLUME(ng) % rhoAm(K)
  ENDDO
  close(unit=15)
END PROGRAM
