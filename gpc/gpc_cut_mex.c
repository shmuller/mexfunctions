/* res = gpc_cut_mex(xy,n)
 * Cut list of polygons using the General Polygon Clipper (GPC) library.
 *
 * URL: 
 * http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 -Igpc232 gpc_cut_mex.c gpc232/gpc.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 -Igpc232 gpc_cut_mex.c gpc232/gpc.c
 *
 * S. H. Muller, 2010/09/04
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "gpc.h"

void mxExport(mxArray **L, gpc_polygon *p, int J)
{
    mwSize dims[2];
    int *m, *M, Mc, *n, *N, Nc, i, j;
    double *d;
    
    dims[0] = J;
    L[2] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    M = (int*) mxGetPr(L[2]);
    for(j=0,m=M,Mc=0; j<J; ++j,++m) {
        Mc += *m = p[j].num_contours;
    }
    dims[0] = Mc;
    L[1] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    N = (int*) mxGetPr(L[1]);
    for(j=0,n=N,Nc=0; j<J; ++j) {
        for(i=0; i<M[j]; ++i,++n) {
            Nc += *n = p[j].contour[i].num_vertices;
        }
    }
    dims[0] = 2; dims[1] = Nc;
    L[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    d = mxGetPr(L[0]);
    for(j=0,n=N; j<J; ++j) {
        for(i=0; i<M[j]; ++i,++n) {
            memcpy(d, p[j].contour[i].vertex, *n*2*sizeof(double));
            d += *n*2;
            if (p[j].hole && p[j].hole[i]) *n = -*n;
        }
    }
}

void mkpoly(gpc_polygon *p, gpc_vertex_list *v)
{
    p->num_contours = 1;
    p->hole = NULL;
    p->contour = v;
}

void cut(gpc_polygon *p, gpc_polygon *q, gpc_polygon *pi)
{
    gpc_polygon pp, qq;
    gpc_polygon_clip(GPC_INT, p, q, pi);
    gpc_polygon_clip(GPC_DIFF, p, q, &pp);
    gpc_polygon_clip(GPC_DIFF, q, p, &qq);
    if (p->hole) gpc_free_polygon(p);
    if (q->hole) gpc_free_polygon(q);
    *p = pp;
    *q = qq;
}

void mexFunction(int nL, mxArray **L, int nR, const mxArray **R)
{
    int i,j;
    const int m = mxGetNumberOfElements(R[1]);
    const int *n, *N = (const int*) mxGetPr(R[1]);
    const double *xy, *XY = mxGetPr(R[0]);
    
    int M = (1 << m) - 1;
    gpc_vertex_list *v, *V = malloc(m*sizeof(gpc_vertex_list));
    gpc_polygon *p, *q, *pi, *P = malloc(M*sizeof(gpc_polygon));
    
    for(i=0,xy=XY,v=V,p=P; i<m; xy+=2*N[i],++v,++p,++i) {
        v->num_vertices = N[i];
        v->vertex = (gpc_vertex*) xy;
    }
    
    for(i=0,v=V,p=P; i<m; ++i,p=pi) {
        mkpoly(p, v++);
        for(q=P,pi=p+1; q<p; ) {
            cut(p, q++, pi++);
        }
    }
        
    printf("%d, %d\n", pi-P, M);
    
    mxExport(L, P, M);
    
    for(i=0,p=P; i<M; ++i,++p) {
        gpc_free_polygon(p);
    }
    free(P);
    free(V);
}
