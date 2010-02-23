/* y = mediansmoothmex(x,w)
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 mediansmoothmex.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 mediansmoothmex.c
 *
 * S. H. Muller, 2010/02/22
 */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

typedef double elem_type;

int compar(const void *a, const void *b)
{
    elem_type d = *(*(const elem_type**)b) - *(*(const elem_type**)a);
    return (d<0)-(d>0);
}

int median(const elem_type **C, const int M, const elem_type *x)
{
    register int i;
    const elem_type **c;
    
    for (i=0,c=C; i<M; i++) {
        *c++ = x+i;
    }
    
    qsort(C,M,sizeof(elem_type*),compar);
    return C[M/2]-x;
}

int median_replace(const elem_type **C, const int M, const elem_type *x, 
    const int del, const int ind)
{
    register int i;
    const elem_type **b, **e;
    elem_type xd = x[del], xi = x[ind];
    
    /* find element to be replaced or position of new element */
    for (i=0,b=C; i<M && *(*b) < xd && *(*b) < xi; i++,b++);
    
    if (*(*b) >= xd) {
        /* found element to be replaced - find new element now */
        for (e=b+1,i++; i<M && *(*e) < xi; i++,e++); e--;
        if (e > b) memmove(b,b+1,(e-b)*sizeof(elem_type*));
        *e = x+ind;
    } else {
        /* found position of new element - find old element now */
        for (e=b; i<M && *(*e) < xd; i++,e++);
        if (e > b) memmove(b+1,b,(e-b)*sizeof(elem_type*));
        *b = x+ind;
    }
    return C[M/2]-x;
}

void print(const elem_type **C, const int M, const elem_type *x)
{
    register int i;
    const elem_type **c;
    
    for (i=0,c=C; i<M; i++,c++) {
        printf("%i: %f\n", *c-x, *(*c));
    }
    printf("\n");
}

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i;
    size_t M;
    
    const mwSize ndims = mxGetNumberOfDimensions(R[0]);
    const mwSize *dims = mxGetDimensions(R[0]);
    const mxClassID mxID = mxGetClassID(R[0]);
    const size_t N = mxGetNumberOfElements(R[0]);
    const size_t bytes = mxGetElementSize(R[0]);
    const double *x = mxGetData(R[0]);
    
    const int *w = mxGetData(R[1]);
    
    L[0] = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int *ind = mxGetData(L[0]);
    ind += *w+1;
    
    M = 2*(*w)+1;
    const elem_type **C = malloc(M*sizeof(elem_type*));
    
    ind[-1] = median(C,M,x);
    print(C,M,x);
    
    for (i=0; i<N-M; i++) {
        ind[i] = median_replace(C,M,x,i,i+M);
        /* ind[i] = i+1 + median(C,M,x+i+1); */
        print(C,M,x);
    }
    
    free(C);
    
    /*
    if (*w+1 >= N) {
        int t = median(x,N);
        for (i=0; i<N; i++) *ind++ = t;
    } else {
        for (i=0; i<*w; i++) {
            *ind++ = median(x, i+1+*w);
        }
        for (i=0; i<N-2*(*w); i++) {
            *ind++ = i + median(x+i*bytes, 2*(*w)+1);
        }
        for (; i<N-*w; i++) {
            *ind++ = i + median(x+i*bytes, N-i);
        }
    }
    */
}
