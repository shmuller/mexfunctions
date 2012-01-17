#include "string.h"

#include "gauss_legendre.h"

#include "quad.h"

void looppar(const int rank, const int *N, int *s, int *no, int *ni, const int dim)
{
    int i, j=dim-1, r=rank-1;
    
    for(i=0,*ni=1;i<j;i++) *ni *= N[i];        // length of inner loop
    for(i=r,*no=1;i>j;i--) *no *= N[i];        // length of outer loop
	
    // stride of outer loop (minus 1 due to start at end of inner loop)
    *s = *ni*(N[j]-1);
}

void filldim(int rank, const int *N, double *X, const double *x, int dim)
{
    int i,j,k,s,no,ni,nj=N[dim-1];
    const double *p;
    double *q, xj;
    
    looppar(rank,N,&s,&no,&ni,dim);
    
    for(k=no,q=X; k--; ) for(j=nj,p=x; j--; ) for(i=ni,xj=*p++; i--; ) {
        *q++ = xj;
    }
}

void filldim2(int rank, const int *N, double *X, const double *x, int dim)
{
    int i,j,k,s,no,ni,nj=N[dim-1],siz;
    double *q, xj;
    
    looppar(rank,N,&s,&no,&ni,dim);
    s = ni*nj;
    siz = s*sizeof(double);
    
    for(j=nj,q=X; j--; ) for(i=ni,xj=*x++; i--; ) {
        *q++ = xj;
    }
    for(k=no-1; k--; q+=s) {
        memcpy(q,X,siz);
    }
}

void weighted_sum(double *y, const double *w, int n, int ni)
{
    int i, k, m=(n+1)>>1, o=n&1;
    double *p, *q, wi;
    
    // add symmetric pairs for outermost dimension
    for(i=m-o,p=y+o*ni,q=y+m*ni; i--; ) {
        for(k=ni; k--; ) {
            *p++ += *q++;
        }
    }
    
    // store weighted sum on first element
    for(k=ni,q=y,wi=*w++; k--; ) {
        *q++ *= wi;
    }
    for(i=m-1; i--; ) {
        for(k=ni,p=y,wi=*w++; k--; ) {
            *p++ += wi*(*q++);
        }
    }
}
    
double *scaleX(int n, double *x, double A, double B)
{
    int i, m=(n+1)>>1, o=n&1;
    double Ax, *p, *q, *X = malloc(n*sizeof(double));
    
    p = X; q = X+m;
    if (o) {
        *p++ = B;
    }
    for(i=o; i<m; i++) {
        Ax = A*x[i];
        *p++ = B+Ax;
        *q++ = B-Ax;
    }
    return X;
}

double getScaledX(int N, int *L, double *X, double *ab, int *dtbl, double **x, double **w, int d)
{
    double A = 0.5*(ab[1]-ab[0]);
	double B = 0.5*(ab[1]+ab[0]);
     
    // load abscissae and weight table
    *dtbl = gauss_legendre_load_tbl(L[d], x, w);
    
    // populate X with scaled abscissae
    double *xx = scaleX(L[d], *x, A, B);
    
    filldim2(N, L, X, xx, d+1);
    
    free(xx);
    
    return A;
}
