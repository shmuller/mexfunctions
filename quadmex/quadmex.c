/* y = quadmex(fun, x)
 * Use gauss_legendre() to calculate integral over Matlab function fun.
 *
 * Compile mex file with (using gnumex):
 *
 * mex -v quadmex.c gauss_legendre.o
 *
 * S. H. Muller, 2012/01/12
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "string.h"

#include "gauss_legendre.h"

typedef struct {
    int M;
    int N;
    int L;
    double *data;
} Data;

void initData(Data *D, int M, int N, double *data)
{
    D->M = M;
    D->N = N;
    D->L = M*N;
    D->data = data;
}

void getData(Data *D, mxArray *R)
{
    D->M = mxGetM(R);
    D->N = mxGetN(R);
    D->L = D->M*D->N;
    D->data = mxGetData(R);
}

void setData(mxArray *R, Data *D)
{
    mxSetM(R, D->M);
    mxSetN(R, D->N);
    mxSetData(R, D->data);
}

void looppar(const int rank, const int *N, int *s, int *no, int *ni, const int dim)
{
    int i, j=dim-1, r=rank-1;
    
    for(i=0,*ni=1;i<j;i++) *ni *= N[i];        // length of inner loop
    for(i=r,*no=1;i>j;i--) *no *= N[i];        // length of outer loop
	
    // stride of outer loop (minus 1 due to start at end of inner loop)
    *s = *ni*(N[j]-1);
}
 
void filldim(const int rank, const int *N, double *X, const double *x, const int dim)
{
    int i,j,k,s,no,ni;
    double *p, *q, xi;
    
    looppar(rank,N,&s,&no,&ni,dim);
    
    for(i=N[dim-1],q=X; i--; q+=ni) {
        xi = *x++;
        for(j=no,p=q; j--; p+=s) for(k=ni; k--; ) {
            *p++ = xi;
        }
    }
}

double *scaleX(int n, double *x, double A, double B)
{
    int i, m=(n+1)>>1, o=n&1;
    double Ax, *p, *q, *X = malloc(n*sizeof(double));
    
    p = X+m; q = p-1;
    if (o) {
        *q-- = B;
    }
    for(i=o; i<m; i++) {
        Ax = A*x[i];
        *p++ = B+Ax;
        *q-- = B-Ax;
    }
    return X;
}

double getScaledX(int n, double *X, double *ab, int *dtbl, double **x, double **w, int M)
{
    int i, j, m=(n+1)>>1, o=n&1, M2=M*2, Mo=M*o, Mm = M*m;
    double A, B, Ax, BpAx, BmAx, *p, *q;
    
    A = 0.5*(ab[1]-ab[0]);
	B = 0.5*(ab[1]+ab[0]);
     
    // load abscissae and weight table
    *dtbl = gauss_legendre_load_tbl(n, x, w);
    
    // populate X with scaled abscissae
    double *xx = scaleX(n, *x, A, B);
    int N[] = {M,n};
    
    filldim(2,N,X,xx,2);
    
    free(xx);
    
    return A;
}


mxArray *gauss_legendre_matlab(int n, int d, int nR, mxArray **R)
{
    mxArray *L, *Y;
    
	double* x = NULL;
	double* w = NULL;
	double A, wi;
	int i, j, k, dtbl, o, m, M, M2, Mo, Mm, D = nR-1;

    o = n&1;
	m = (n+1)>>1;

    double *p, *q, *y;
    
    Data *ND = malloc(2*D*sizeof(Data));
    Data *OD = ND+D, *od, *nd;
    
    int Lm, Ln, Lmn;
    
    // find total lengths of parameter and integral dimensions
    for(i=1,od=OD,Lm=Ln=1; i<=D; i++,od++) {
        getData(od, R[i]);
        if (i == d) {
            Ln *= n;
        } else {
            Lm *= od->L;
        }
    }
    Lmn = Lm*Ln;
    
    int N[] = {Lm,Ln};
    
    double *mem = malloc(Lmn*D*sizeof(double));
    
    // generate tensor product arguments
    for(i=1,od=OD,nd=ND,p=mem; i<=D; i++,od++,nd++,p+=Lmn) {
        if (od->L == 1) {
            continue;  // do not expand scalar parameter
        }
        initData(nd, Lm, Ln, p);
        if (i == d) {
            A = getScaledX(n, nd->data, od->data, &dtbl, &x, &w, Lm);
        } else {
            filldim(2,N,nd->data,od->data,1);
        }  
        setData(R[i], nd);
    }
    
    // evaluate function at tensor product arguments
    mexCallMATLAB(1, &Y, nR, R, "feval");
    y = mxGetData(Y);
    
    // re-attach old data
    for(i=1,od=OD; i<=D; i++,od++) {
        setData(R[i], od);   
    }
    
    int s, no, ni;
    double *P, *Q, *W;
    
    looppar(2,N,&s,&no,&ni,2);
    
    // add symmetric pairs
    for(i=m-o,P=y,Q=y+(n-1)*ni; i--; P+=ni,Q-=ni) {
        for(j=no,p=P,q=Q; j--; p+=s,q+=s) for(k=ni; k--; ) {
            *q++ += *p++;
        }
    }
    
    // create output matrix
    M = mxGetM(Y);
    L = mxCreateNumericMatrix(M,1,mxDOUBLE_CLASS,mxREAL);
    P = mxGetData(L);
    
    // populate output matrix with weighted sums
    for(i=m,Q=y+(m-o)*ni,W=w; i--; Q+=ni,W++) {
        for(j=no,p=P,q=Q,wi=*W; j--; q+=s) for(k=ni; k--; ) {
            *p++ += wi*(*q++);
        }
    }
    
    // apply final scaling
    for(j=0,p=P; j<M; j++) {
        *p++ *= A;
    }
    
    mxDestroyArray(Y);
    free(mem);
    free(ND);
    
	if (dtbl) {
		free(x);
		free(w);
	}
	return L;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const mwSize ndims = mxGetNumberOfDimensions(R[1]);
    const mwSize *dims = mxGetDimensions(R[1]);
    const int npts = mxGetNumberOfElements(R[1]);
    
    int n = *((int*)mxGetData(R[0]));
    int d = *((int*)mxGetData(R[1]));
    
    L[0] = gauss_legendre_matlab(n, d, nR-2, (mxArray**)(R+2));
}

