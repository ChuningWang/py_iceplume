MODULE mod_iceplume_py
  USE mod_kinds
  implicit none
  integer, parameter :: Ngrids = 1
  integer, parameter :: itemp = 1
  integer, parameter :: isalt = 2
!
  integer, parameter :: Nsrc(1) = [1]
  integer, parameter :: NT(1) = [3]
  integer, parameter :: N(1) = [40]
  real(r8), parameter :: dt(1) = [30.0d0]
!
  logical, parameter :: usePotTemp = .true.
END MODULE
