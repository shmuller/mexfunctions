/* y = quadmex(int32([di,dj]),int32([ni,nj]),@fun,p1,...,li,...,lj,...,pn)
 * Use gauss_legendre() to calculate integral over Matlab function fun.
 *
 * Compile mex file with (using gnumex):
 *
 * mex -v quadmex.c quad.o gauss_legendre.o
 *
 * S. H. Muller, 2012/01/12
 */

#include "mex.h"
#include "math.h"
#include "string.h"

#include "quad.h"

#define max(a,b) ((b) > (a) ? (b) : (a))

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

mxArray *gauss_legendre_matlab(int nI, const int *d, const int *n, int nR, mxArray **R)
{
    mxArray *res, *Y;
    
	double* x = NULL;
	double* w = NULL;
	double A;
	int i, c, dtbl, D = nR-1;

    const int *dd;
    double *mem, *p, *q, *y;
    
    Data *ND = malloc(2*D*sizeof(Data));
    Data *OD = ND+D, *od, *nd;
    
    int *L = malloc((nI+1)*sizeof(int));
    int Lm, Ln, Lmn, ni;
    
    // find total lengths of parameter and integral dimensions
    for(i=1,od=OD,Lm=Ln=1,c=0; i<=D; i++,od++) {
        getData(od, R[i]);
        if (c < nI && i == d[c]) {
            L[c+1] = n[c];
            Ln *= n[c];
            ++c;
        } else {
            Lm = max(Lm,od->L);
        }
    }
    L[0] = Lm;
    Lmn = Lm*Ln;
        
    if ((mem=malloc(Lmn*D*sizeof(double)))==NULL) {
        free(L);
        free(ND);
        mexErrMsgTxt("Out of memory");
    }
    
    // generate tensor product arguments
    for(i=1,od=OD,nd=ND,p=mem,A=1.0,c=0; i<=D; i++,od++,nd++,p+=Lmn) {
        if (od->L == 1) {
            continue;  // do not expand scalar parameter
        }
        initData(nd, Lm, Ln, p);
        if (c < nI && i == d[c]) {
            A *= getScaledX(nI+1, L, nd->data, od->data, &dtbl, &x, &w, c+1);
            ++c;
        } else {
            filldim2(nI+1, L, nd->data, od->data, 1);
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
    
    // calculated weighted sums
    for(i=c,ni=Lmn; i>0; i--) {
        ni/=L[i];
        weighted_sum(y,w,L[i],ni);
    }

    // create output matrix
    res = mxCreateNumericMatrix(Lm,1,mxDOUBLE_CLASS,mxREAL);
    
    // apply final scaling
    p = mxGetData(res);
    for(i=Lm,q=y; i--; ) {
        *p++ = A*(*q++);
    }
    
    mxDestroyArray(Y);
    free(mem);
    free(L);
    free(ND);
    
	if (dtbl) {
		free(x);
		free(w);
	}
	return res;
}


void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    const int nI = mxGetNumberOfElements(R[0]);
    const int *d = mxGetData(R[0]);
    const int *n = mxGetData(R[1]);
    
    L[0] = gauss_legendre_matlab(nI, d, n, nR-2, (mxArray**)(R+2));
}
