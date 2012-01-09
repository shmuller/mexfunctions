% How to compile on Windows with standalone lapack:
%
% (1) Build lapack and blas on MSYS with MINGW gfortran.
% (2) Run gnumex.
% (3) On D:\matlab\gnumex\linkmex.pl, change "$linker = 'gcc -shared';" to
%
%     $linker = 'gfortran -shared';
%
% (4) Run make.m.

%mex -v OPTIMFLAGS=-O3 -Ddgesv=dgesv_ triblockmex.c "C:\MinGW32-xy\lib\liblapack.a" "C:\MinGW32-xy\lib\libblas.a"

mex -v OPTIMFLAGS=-O3 triblockmex.c -L"C:\PROGRA~1\MATLAB\R2007a\extern\lib\win32\lcc" -lmwlapack

%mex -v OPTIMFLAGS=-O3 -Ddgetrf=dgetrf_ -Ddgetrs=dgetrs_ -Ddgemm=dgemm_ triblockmex.c "lib\liblapack.a" "lib\libcblas.a" "lib\libf77blas.a" "lib\libatlas.a"