MODULE mod_py_iceplume
  USE mod_kinds
  implicit none
  integer, parameter :: Ngrids = 1
  integer, parameter :: itemp = 1
  integer, parameter :: isalt = 2
!
  integer, parameter :: Nsrc(1) = [1]
  integer, parameter :: NT(1) = [3]
  integer            :: N(1)
  real(r8)           :: dt(1)
!
  logical, parameter :: usePotTemp = .true.
END MODULE
