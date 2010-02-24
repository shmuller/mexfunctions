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

int compar(const void *a, const void *b)
{
    double d = *(*(const elem_ptr*)b) - *(*(const elem_ptr*)a);
    return (d<0)-(d>0);
}

int median(elem_ptr *C, const int M, const double *x)
{
    register int i;
    elem_ptr *c;
    
    for (i=0,c=C; i<M; i++,c++) {
        *c = x+i;
    }
    
    qsort(C,M,sizeof(elem_ptr),compar);
    return C[M/2]-x;
}

int median_replace(elem_ptr *C, const int M, const elem_ptr x, 
    const int del, const int ind)
{
    register int i;
    elem_ptr *b, *e, *pd, *pi, *pb, *pe;
    elem_type xi = x[ind], xd = x[del], xm = (xi < xd) ? xi : xd;
    
    /*
    for (i=0,pd=C; i<M && *(*pd) < xd; i++,pd++);
    for (i=0,pi=C; i<M && *(*pi) < xi; i++,pi++);
    */
    for (pd=C; *(*pd) < xd; pd++);
    for (pi=C; *(*pi) < xi; pi++);
    
    if (pd < pi) {
        pi--;
        pb = pd; pe = pi;
    } else {
        pb = pi; pe = pd;
    }
    
    /* find element to be replaced or position of new element */
    for (i=0,b=C; i<M && *(*b) < xm; i++,b++);
    
    if (*(*b) == xd) {
        /* found element to be replaced - find new element now */
        for (e=b+1,i++; i<M && *(*e) < xi; i++,e++); e--;
        if (e > b) memmove(b,b+1,(e-b)*sizeof(elem_ptr));
        *e = x+ind;
    } else {
        /* found position of new element - find old element now */
        for (e=b; i<M && *(*e) < xd; i++,e++);
        if (e > b) memmove(b+1,b,(e-b)*sizeof(elem_ptr));
        *b = x+ind;
    }
    
    if (pb != b || pe != e) {
        mexErrMsgTxt("Error calculating memmove block!!!");
    }
    
    return C[M/2]-x;
}

void print(const elem_ptr *C, const int M, const elem_ptr x)
{
    register int i;
    const elem_ptr *c;
    
    for (i=0, c=C; i<M; i++,c++) {
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
    
    M = 2*(*w)+1;
    elem_type sent[] = {-111., 111.};
    elem_ptr *C = (elem_ptr*) malloc((M+2)*sizeof(elem_ptr)) + 1;
    C[-1] = &sent[0]; C[M] = &sent[1];
    
    ind += *w;
    *ind++ = median(C,M,x);
    print(C,M,x);
    
    for (i=0; i<N-M; i++) {
        *ind++ = median_replace(C,M,x,i,i+M);
        /* *ind++ = i+1 + median(C,M,x+i+1); */
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
