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

	const NDArray AA = args(0).array_value();
	const NDArray BB = args(1).array_value();
	const int32NDArray IP = args(2).int32_array_value();
	const NDArray CC = args(3).array_value();
	const NDArray DD = args(4).array_value();

	double *A = (double*) AA.data();
	double *B = (double*) BB.data();
	int *IPIV = (int*) IP.data();
	double *C = (double*) CC.data();
	double *D = (double*) DD.data();

	const dim_vector dv = BB.dims();
	int N = dv(0), NRHS = dv(1);

	int LDA = N, LDB = N, INFO = 0, N2 = N/2, NN2 = N*N2;
	char TRANSN = 'N', TRANST = 'T';
	double ALPHA = -1.0, BETA = 1.0;

	dgetrf(&N, &N, A, &LDA, IPIV, &INFO);
	
	if (INFO == 0) {
 		dgetrs(&TRANST, &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
 	}
        
 	D += NN2;
 	dgemm(&TRANST, &TRANSN, &N, &N2, &N, &ALPHA, B, &N, C, &N, &BETA, D, &N);
	
	return retval;
}


