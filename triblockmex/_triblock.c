#include "_triblock.h"

extern void dgesv(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
    double *B, int *LDB, int *INFO);

extern void dgetrf(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

extern void dgetrs(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
    int *IPIV, double *B, int *LDB, int *INFO);

extern void dgemm(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
    double *ALPHA, double *A, int *LDA, double *B, int *LDB,
    double *BETA, double *C, int *LDC);


void _triblock(double *D1, double *U2, int *IPIV, double *L1, double *D2, int is_half,
    int N, int NRHS, int *INFO)
{
    char TRANSN = 'N', TRANST = 'T';
	double ALPHA = -1.0, BETA = 1.0;

	dgetrf(&N, &N, D1, &N, IPIV, INFO);
	
    if (!is_half) 
    {
        // U2 = D1\U2 
        // D2 = D2 - L1*U2
        dgetrs(&TRANSN, &N, &NRHS, D1, &N, IPIV, U2, &N, INFO);
        
        dgemm(&TRANSN, &TRANSN, &N, &NRHS, &N, &ALPHA, L1, &N, U2, &N, &BETA, D2, &N);
    } 
    else 
    {
        dgetrs(&TRANST, &N, &NRHS, D1, &N, IPIV, U2, &N, INFO);
        
        int N2 = N/2, NN2 = N*N2;
 	    D2 += NN2;
 	    dgemm(&TRANST, &TRANSN, &N, &N2, &N, &ALPHA, U2, &N, L1, &N, &BETA, D2, &N);
    }
}

