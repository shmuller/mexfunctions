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

typedef struct {
    int ind;
    const double *elem;
} container;

int compar(const void *a, const void *b)
{
    double d = *((const container*)b)->elem - *((const container*)a)->elem;
    return (d<0)-(d>0);
}

void replace(container *C, const int M, const double *x, 
    const int del, const int ind)
{
    register int i;
    container *b, *e;
    double xi = x[ind];
    
    /* find element to be replaced or position of new element */
    for (i=0,b=C; i<M && b->ind != del && *b->elem < xi; i++,b++);
    
    if (b->ind == del) {
        /* found element to be replaced - find new element now */
        for (e=b+1,i++; i<M && *e->elem < xi; i++,e++); e--;
        if (e > b) memmove(b,b+1,(e-b)*sizeof(container));
        e->ind = ind;
        e->elem = x+ind;
    } else {
        /* found position of new element - find old element now */
        for (e=b; i<M && e->ind != del; i++,e++);
        if (e > b) memmove(b+1,b,(e-b)*sizeof(container));
        b->ind = ind;
        b->elem = x+ind;
    }
}

void print(const container *C, const int M)
{
    register int i;
    const container *c;
    
    for (i=0, c=C; i<M; i++,c++) {
        printf("%i: %f\n", c->ind, *c->elem);
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
    
    M = 2*(*w)+1;
    container *C, *c;
    C = malloc(M*sizeof(container));
    
    for (i=0, c=C; i<M; i++,c++) {
        c->ind = i;
        c->elem = x+i;
    }
    
    ind += *w+1;
    
    qsort(C,M,sizeof(container),compar);
    ind[-1] = C[M/2].ind;
    
    /* print(C,M); */
    
    for (i=0; i<N-M; i++) {
        replace(C,M,x,i,i+M);
        ind[i] = C[*w].ind;
        
        /* print(C,M); */
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
