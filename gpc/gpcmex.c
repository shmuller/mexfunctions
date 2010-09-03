/* res = gpcmex(x1,y1,x2,y2)
 * MEX interface for the General Polygon Clipper (GPC) library.
 *
 * URL: 
 * http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
 *
 * Compile on Linux: 
 * mex -v COPTIMFLAGS=-O3 LDOPTIMFLAGS=-O3 -Igpc232 gpcmex.c gpc232/gpc.c
 * 
 * Compile on Windows
 * mex -v OPTIMFLAGS=-O3 -Igpc232 gpcmex.c gpc232/gpc.c
 *
 * S. H. Muller, 2010/09/02
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "gpc.h"

void mxExport(mxArray **L, gpc_polygon *p)
{
    mwSize dims[2];
    double *d;
    int *n, nc, N, i;
    
    dims[0] = N = p->num_contours;
    L[1] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    n = (int*) mxGetPr(L[1]);
    for(i=0,nc=0; i<N; i++) {
        nc += n[i] = p->contour[i].num_vertices;
    }
    dims[0] = 2; dims[1] = nc;
    L[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    for(i=0,d=mxGetPr(L[0]); i<N; d+=2*n[i],i++) {
        memcpy(d, p->contour[i].vertex, 2*n[i]*sizeof(double));
    }
}

void mexFunction(int nL, mxArray **L, int nR, const mxArray **R)
{
    const int n1 = mxGetNumberOfElements(R[0])/2;
    const int n2 = mxGetNumberOfElements(R[1])/2;
    
    gpc_vertex_list l1 = {n1, (gpc_vertex*) mxGetPr(R[0])};
    gpc_vertex_list l2 = {n2, (gpc_vertex*) mxGetPr(R[1])};
    
    gpc_polygon p1 = {1, NULL, &l1}, p2 = {1, NULL, &l2}, p;
    
    gpc_polygon_clip(GPC_INT, &p1, &p2, &p);
    mxExport(L, &p);
    gpc_free_polygon(&p);
    
    if (nL == 6) {
        gpc_polygon_clip(GPC_DIFF, &p1, &p2, &p);
        mxExport(L+2, &p);
        gpc_free_polygon(&p);
    
        gpc_polygon_clip(GPC_DIFF, &p2, &p1, &p);
        mxExport(L+4, &p);
        gpc_free_polygon(&p);
    }    
}
