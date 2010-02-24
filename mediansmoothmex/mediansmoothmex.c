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
        printf("%03d: %f\n", *c-x, *(*c));
    }
    printf("\n");
}

int compar(const void *a, const void *b)
{
    elem_type d = *(*(const elem_ptr*)b) - *(*(const elem_ptr*)a);
    return (d<0)-(d>0);
}

int median_init(elem_ptr * const C, const int M, const elem_ptr x)
{
    register int i;
    elem_ptr *c;
    
    for (i=0,c=C; i<M; i++) {
        *c++ = x+i;
    }
    
    qsort(C,M,sizeof(elem_ptr),compar);
    return C[(M-1)/2]-x;
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

int median_add(elem_ptr * const C, const int M, const elem_ptr x, 
    const int ind, elem_ptr **save)
{
    elem_ptr *pi = save[1];
    
    pi = find_spot(pi,x[ind]);
    memmove(pi+1,pi,(M-(pi-C))*sizeof(elem_ptr));
    *pi = x+ind;
    
    save[1] = pi;
    return C[(M-1)/2]-x;
}

int median_remove(elem_ptr * const C, const int M, const elem_ptr x, 
    const int del, elem_ptr **save)
{
    elem_ptr *pd = save[0];
    
    pd = find_spot(pd,x[del]);
    memmove(pd,pd+1,(M-(pd-1-C))*sizeof(elem_ptr));
    
    save[0] = pd;
    return C[(M-1)/2]-x;
}

int median_replace(elem_ptr * const C, const int M, const elem_ptr x, 
    const int del, const int ind, elem_ptr **save)
{
    elem_ptr *pd = save[0], *pi = save[1];
    
    pd = find_spot(pd,x[del]);
    pi = find_spot(pi,x[ind]);
    
    if (pd < pi) {
        pi--;
        if (pi > pd) memmove(pd,pd+1,(pi-pd)*sizeof(elem_ptr));
    } else if (pi < pd) {
        memmove(pi+1,pi,(pd-pi)*sizeof(elem_ptr));
    }
    
    if (pi < C || pi >= C+M) {
        printf("Error: Index out of bounds!\n");
        exit(1);
    }
    
    *pi = x+ind;
    
    //save[0] = pd;
    //save[1] = pi;
    
    return C[(M-1)/2]-x;
}

int median(const elem_ptr x, const size_t N, const int w, int *ind)
{
    register int i;
    const size_t M = 2*w+1;
    
    //elem_type sent[] = {-DBL_MAX, DBL_MAX};
    elem_type sent[] = {-111., 111.};
    elem_ptr * const C = malloc((M+2)*sizeof(elem_ptr)) + 1;
    C[-1] = &sent[0]; C[w+1] = C[M] = &sent[1];
    
    elem_ptr *save[] = {C,C};
    /*
    printf("Initialize...\n");
    *ind++ = median_init(C,w+1,x);
    print(C,w+1,x);
    
    printf("Add...\n");
    for (i=0; i<w; i++) {
        *ind++ = median_add(C,w+2+i,x,w+1+i,save);
        print(C,w+2+i,x);
    }
    */
    ind += w;
    *ind++ = median_init(C,M,x);
    print(C,M,x);
    
    //printf("Replace...\n");
    for (i=0; i<N-M; i++) {
        *ind++ = median_replace(C,M,x,i,i+M,save);
        print(C,M,x);
    }
    
    /*
    printf("Remove...\n");
    for (i=0; i<w; i++) {
        *ind++ = median_remove(C,M-i-1,x,N-M+i,save);
        print(C,M-i-1,x);
    }
    */
    free(C-1);
    return 0;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const mwSize ndims = mxGetNumberOfDimensions(R[0]);
    const mwSize *dims = mxGetDimensions(R[0]);
    const mxClassID mxID = mxGetClassID(R[0]);
    const size_t N = mxGetNumberOfElements(R[0]);
    const size_t bytes = mxGetElementSize(R[0]);
    const double *x = mxGetData(R[0]);
    
    const int *w = mxGetData(R[1]);
    
    L[0] = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int *ind = mxGetData(L[0]);
    
    median(x,N,*w,ind);
}
