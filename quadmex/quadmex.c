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
    double *p, *q, xk;
    
    looppar(rank,N,&s,&no,&ni,dim);
    
    for(k=N[dim-1],q=X; k--; q+=ni) {
        xk = *x++;
        for(j=no,p=q; j--; p+=s) for(i=ni; i--; ) {
            *p++ = xk;
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
    
    /*
    if(o) {
        q = X+Mm-M;
        for(j=0; j<M; j++) {
            *q++ = B;
        }
    }
    for(i=o,p=X+Mm,q=p-M-Mo; i<m; i++,q-=M2) {
        Ax = A*(*x)[i];
        BpAx = B+Ax;
        BmAx = B-Ax;
        for(j=0; j<M; j++) {
            *p++ = BpAx;
            *q++ = BmAx;
        }
    }
     */
    
    return A;
}


mxArray *gauss_legendre_matlab(int n, int d, int nR, mxArray **R)
{
    mxArray *L, *Y;
    
	double* x = NULL;
	double* w = NULL;
	double A, wi;
	int i, j, dtbl, o, m, M, M2, Mo, Mm, D = nR-1;

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
    double *mem = malloc(Lmn*D*sizeof(double));
    
    //printf("D = %d, Lm = %d, Ln = %d\n", D, Lm, Ln);
    
    // generate tensor product arguments
    for(i=1,od=OD,nd=ND,p=mem; i<=D; i++,od++,nd++,p+=Lmn) {
        if (od->L == 1) {
            continue;  // do not expand scalar parameter
        }
        initData(nd, Lm, Ln, p);
        if (i == d) {
            A = getScaledX(n, nd->data, od->data, &dtbl, &x, &w, Lm);
            //printf("A = %f\n", A);
        } else {
            for(j=0,q=nd->data; j<Ln; j++,q+=Lm) {
                memcpy(q, od->data, Lm*sizeof(double));
            }
        }
        /*
        for(j=0; j<Lmn; j++) {
            printf("%f\n", nd->data[j]);
        }
        printf("\n");
        */  
        setData(R[i], nd);
    }
    
    //return R[d];
    
    // evaluate function at tensor product arguments
    mexCallMATLAB(1, &Y, nR, R, "feval");
    
    // re-attach old data
    for(i=1,od=OD; i<=D; i++,od++) {
        setData(R[i], od);   
    }
    
    M = mxGetM(Y);
    y = mxGetData(Y);
    
    M2 = M*2;
    Mo = M*o;
    Mm = M*m;
    
    // add symmetric pairs
    for(i=o,p=y+Mm,q=p-M-Mo; i<m; i++,q-=M2) {
        for(j=0; j<M; j++) {
            *p++ += *q++;
        }
    }
    
    // create output matrix
    L = mxCreateNumericMatrix(M,1,mxDOUBLE_CLASS,mxREAL);
    p = mxGetData(L);
    
    // populate output matrix with weighted sums
    for(i=0,q=y+Mm-Mo; i<m; i++,p-=M) {
        for(j=0,wi=w[i]; j<M; j++) {
            *p++ += wi*(*q++);
        }
    }
    // apply final scaling
    for(j=0; j<M; j++) {
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

