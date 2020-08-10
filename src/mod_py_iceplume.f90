MODULE mod_py_iceplume
! =====================================================================!
!                                                                      !
! These are the module functions of py-iceplume along.                 !
!                                                                      !
! =====================================================================!
!
! This module stores all global variables.
!
  USE mod_kinds
  implicit none
  integer, parameter :: Ngrids  = 1
  integer, parameter :: itemp   = 1
  integer, parameter :: isalt   = 2
!
  integer, parameter :: Nsrc(1) = [3]
  integer, parameter :: NT(1)   = [3]
  integer            :: N(1)
  real(r8)           :: dt
!
  logical, parameter :: usePotTemp = .true.
END MODULE
