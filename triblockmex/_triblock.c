#include "_triblock.h"

extern void dgesv(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
    double *B, int *LDB, int *INFO);

extern void dgetrf(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

extern void dgetrs(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
    int *IPIV, double *B, int *LDB, int *INFO);

extern void dtrsm(char *SIDE, char *UPLO, char* TRANSA, char *DIAG,
    int *M, int *N, double *ALPHA, double *A, int *LDA, double *B, int *LDB);

extern void dlaswp(int *N, double *A, int *LDA, int *K1, int *K2, int *IPIV, int *INCX);

extern void dgemm(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
    double *ALPHA, double *A, int *LDA, double *B, int *LDB,
    double *BETA, double *C, int *LDC);


void swap_columns(int *IPIV, int n, double *x, int m)
{
    int i, j, k;
    double tmp, *a, *b;

    --x;
    a = x+n*m;
    for(j=n; j--; ) {
        if ((k=IPIV[j]) != j+1) {
            b = x+k*m;
            for(i=m; i--; ) {
                tmp = *a; *a-- = *b; *b-- = tmp;
            }
        } else {
            a -= m;
        }
    }
}


void _triblock(double *D1, double *U2, int *IPIV, double *L1, double *D2, int is_half,
    int N, int NRHS, int *INFO)
{
    int tmp;
    char TRANSN='N', TRANST='T';
    char LEFT='L', RIGHT='R', UPPER='U', LOWER='L', UNIT='U', NONUNIT='N';
    double ALPHA = -1.0, BETA = 1.0, ONE = 1.0;
    int IONE = 1, IMINUSONE = -1;

    if (is_half) {
        tmp = N; N = NRHS; NRHS = tmp;
    }

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
        // U2 = U2/D1
        // D2 = D2 - U2*L1
        dtrsm(&RIGHT, &UPPER, &TRANSN, &NONUNIT, &NRHS, &N, &ONE, D1, &N, U2, &NRHS);

        dtrsm(&RIGHT, &LOWER, &TRANSN, &UNIT, &NRHS, &N, &ONE, D1, &N, U2, &NRHS);

        swap_columns(IPIV, N, U2, NRHS);

        int N2 = N/2, NN2 = N*N2;
        D2 += NN2;
        dgemm(&TRANSN, &TRANSN, &N, &N2, &N, &ALPHA, U2, &N, L1, &N, &BETA, D2, &N);
    }
    
    /*
    else
    {        
        dtrsm(&LEFT, &UPPER, &TRANST, &NONUNIT, &N, &NRHS, &ONE, D1, &N, U2, &N);

        dtrsm(&LEFT, &LOWER, &TRANST, &UNIT, &N, &NRHS, &ONE, D1, &N, U2, &N);

        dlaswp(&NRHS, U2, &N, &IONE, &N, IPIV, &IMINUSONE);

        //dgetrs(&TRANST, &N, &NRHS, D1, &N, IPIV, U2, &N, INFO);
        
        //int N2 = N/2, NN2 = N*N2;
        //D2 += NN2;
        //dgemm(&TRANST, &TRANSN, &N, &N2, &N, &ALPHA, U2, &N, L1, &N, &BETA, D2, &N);
    }
    */

}

