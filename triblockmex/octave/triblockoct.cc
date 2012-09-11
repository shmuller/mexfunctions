/* triblockmex(A,B,IPIV,C,D)
 * Tri-block-diagnonal solve step. Compile with:
 *
 * mkoctfile triblockoct.cc _triblock.c -llapack -lblas
 *
 * ln -s triblockoct.oct triblockmex.oct
 *
 * S. H. Muller, 2012/02/15
 */

#include <octave/oct.h>

#include "../_triblock.h"

DEFUN_DLD(triblockmex, args, nargout, "Tri-block-diagonal solve step")
{
	octave_value retval;

	const NDArray arr_D1 = args(0).array_value();
	const NDArray arr_U2 = args(1).array_value();
	const int32NDArray arr_IPIV = args(2).int32_array_value();
	const NDArray arr_L1 = args(3).array_value();
	const NDArray arr_D2 = args(4).array_value();

	double *D1 = (double*) arr_D1.data();
	double *U2 = (double*) arr_U2.data();
	int *IPIV  = (int*) arr_IPIV.data();
	double *L1 = (double*) arr_L1.data();
	double *D2 = (double*) arr_D2.data();

    int is_half = (args.length() > 5);

	const dim_vector dv = arr_U2.dims();
	int N = dv(0), NRHS = dv(1);
    
	int INFO = 0;
	
    _triblock(D1, U2, IPIV, L1, D2, is_half, N, NRHS, &INFO);

	return retval;
}


