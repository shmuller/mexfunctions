/* X = triblockmex(A,B)
 *
 * S. H. Muller, 2011/12/23
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

extern void dgesv(int *N, int *NRHS, double *A, int *LDA, int *IPIV, 
    double *B, int *LDB, int *INFO);

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{    
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    double *A=mxGetData(R[0]), *B=mxGetData(R[1]);
    int *IPIV=mxGetData(R[2]);
    
    int N = dims[0], NRHS = dims[1], LDA = N, LDB = N, INFO = 0;
    //int *IPIV = malloc(N*sizeof(int));

    //double *a = malloc(N*N*sizeof(double));
    //double *b;
    
    //L[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    //b = mxGetPr(L[0]);
    
    //memcpy(a,A,N*N*sizeof(double));
    //memcpy(b,B,npts*sizeof(double));
    
    dgesv(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
    
    //printf("ndims = %d, N = %d, NRHS = %d, npts = %d, INFO = %d\n", ndims, N, NRHS, npts, INFO);
    
    //free(a);
    //free(IPIV);
    
}
