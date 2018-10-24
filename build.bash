rm *.o *.mod *.exe

gfortran -c ./src/mod_kinds.f90
gfortran -c ./src/mod_iceplume_py.f90
gfortran -c ./src/mod_iceplume.f90
gfortran -c ./src/iceplume_opkd.f90
gfortran -c ./src/iceplume_core.f90
gfortran -c ./src/iceplume_calc.f90
gfortran -c ./src/iceplume.f90

gfortran *.o -o ./iceplume_test.exe

rm *.o *.mod
