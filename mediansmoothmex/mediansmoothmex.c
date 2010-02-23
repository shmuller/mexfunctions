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

typedef const double elem_type;
typedef elem_type* elem_ptr;

void print(const elem_ptr * const C, const int M, const elem_ptr x)
{
    register int i;
    const elem_ptr *c;
    
    for (i=0,c=C; i<M; i++,c++) {
        printf("%i: %f\n", *c-x, *(*c));
    }
    printf("\n");
}

int compar(const void *a, const void *b)
{
    elem_type d = *(*(const elem_ptr*)b) - *(*(const elem_ptr*)a);
    return (d<0)-(d>0);
}

int median(elem_ptr * const C, const int M, const elem_ptr x)
{
    register int i;
    elem_ptr *c;
    
    for (i=0,c=C; i<M; i++) {
        *c++ = x+i;
    }
    
    qsort(C,M,sizeof(elem_ptr),compar);
    return C[M/2]-x;
}

elem_ptr* find_spot(elem_ptr *seed, const elem_type xi)
{
    if (*(*seed) < xi) {
        for (; *(*seed) < xi; seed++);
    } else {
        for (; *(*seed) > xi; seed--);
        if (*(*seed) < xi) seed++;
    }
    return seed;
}

int median_replace(elem_ptr * const C, const int M, const elem_ptr x, 
    const int del, const int ind, elem_ptr **save)
{
    elem_ptr *pd = save[0], *pi = save[1];
    
    pd = find_spot(pd,x[del]);
    pi = find_spot(pi,x[ind]);
    
    if (pd < pi) {
        memmove(pd,pd+1,(pi-pd)*sizeof(elem_ptr));
        pi--;
    } else if (pi < pd) {
        memmove(pi+1,pi,(pd-pi)*sizeof(elem_ptr));
    }
    *pi = x+ind;
    
    save[0] = pd;
    save[1] = pi;
    
    return C[M/2]-x;
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
    elem_type sent[] = {-DBL_MAX, DBL_MAX};
    elem_ptr * const C = malloc((M+2)*sizeof(elem_ptr)) + 1;
    C[-1] = &sent[0]; C[M] = &sent[1];
    
    elem_ptr *save[] = {C,C};
    
    ind[-1] = median(C,M,x);
    print(C,M,x);
    
    for (i=0; i<N-M; i++) {
        ind[i] = median_replace(C,M,x,i,i+M,save);
        /* ind[i] = i+1 + median(C,M,x+i+1); */
        print(C,M,x);
    }
    
    free(C-1);
    
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
