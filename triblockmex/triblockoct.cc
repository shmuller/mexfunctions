/* triblockmex(A,B,IPIV,C,D)
 * Tri-block-diagnonal solve step. Compile with:
 *
 * mkoctfile triblockoct.cc -llapack
 *
 * ln -s triblockoct.oct triblockmex.oct
 *
 * S. H. Muller, 2012/02/15
 */

#include <octave/oct.h>

extern "C"
{
extern void dgesv(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
    double *B, int *LDB, int *INFO);

extern void dgetrf(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

extern void dgetrs(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
    int *IPIV, double *B, int *LDB, int *INFO);

extern void dgemm(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
    double *ALPHA, double *A, int *LDA, double *B, int *LDB,
    double *BETA, double *C, int *LDC);
}


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
	char TRANSN = 'N', TRANST = 'T';
	double ALPHA = -1.0, BETA = 1.0;

	dgetrf(&N, &N, D1, &N, IPIV, &INFO);
	
    if (!is_half) 
    {
        // U2 = D1\U2 
        // D2 = D2 - L1*U2
        dgetrs(&TRANSN, &N, &NRHS, D1, &N, IPIV, U2, &N, &INFO);
        
        dgemm(&TRANSN, &TRANSN, &N, &NRHS, &N, &ALPHA, L1, &N, U2, &N, &BETA, D2, &N);
    } 
    else 
    {
        dgetrs(&TRANST, &N, &NRHS, D1, &N, IPIV, U2, &N, &INFO);
        
        int N2 = N/2, NN2 = N*N2;
 	    D2 += NN2;
 	    dgemm(&TRANST, &TRANSN, &N, &N2, &N, &ALPHA, U2, &N, L1, &N, &BETA, D2, &N);
    }

	return retval;
}


