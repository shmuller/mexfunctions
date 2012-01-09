gcc -c dgesv_ex.c -Ddgesv=dgesv_

#gfortran dgesv_ex.o -o dgesv_ex.exe -llapack -lblas

g77 dgesv_ex.o -o dgesv_ex.exe ./lib/liblapack.a ./lib/libcblas.a ./lib/libf77blas.a ./lib/libatlas.a
