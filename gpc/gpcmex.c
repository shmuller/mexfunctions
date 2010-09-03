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

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i;
    
    const int n1 = mxGetNumberOfElements(R[0])/2;
    const int n2 = mxGetNumberOfElements(R[1])/2;
    
    gpc_vertex_list l1 = {n1, (gpc_vertex*) mxGetPr(R[0])};
    gpc_vertex_list l2 = {n2, (gpc_vertex*) mxGetPr(R[1])};
    
    gpc_polygon p1 = {1, NULL, &l1}, p2 = {1, NULL, &l2}, p3;
    
    int N, *n, nc;
    double *d;
    mwSize dims[2];
    
    gpc_polygon_clip(GPC_INT, &p1, &p2, &p3);
    N = p3.num_contours;
    
    dims[0] = N;
    L[1] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    n = (int*) mxGetPr(L[1]);
    for(i=0,nc=0; i<N; i++) {
        nc += n[i] = p3.contour[i].num_vertices;
    }
    
    dims[0] = 2; dims[1] = nc;
    L[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    for(i=0,d=mxGetPr(L[0]); i<N; d+=2*n[i],i++) {
        memcpy(d, p3.contour[i].vertex, 2*n[i]*sizeof(double));
    }
    
    gpc_free_polygon(&p3);
}
