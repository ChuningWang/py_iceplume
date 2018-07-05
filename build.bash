rm *.o *.mod *.exe

gfortran -c mod_kinds.f90
gfortran -c mod_iceplume_py.f90
gfortran -c mod_iceplume.f90
gfortran -c iceplume_opkd.f90
gfortran -c iceplume_core.f90
gfortran -c iceplume_calc.f90
gfortran -c iceplume.f90

gfortran *.o -o iceplume_test.exe

rm *.o *.mod
