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

gpc_vertex *mk_gpc_vertex(int n, const double *x, const double *y)
{
    register int i;
    gpc_vertex *v = malloc(n*sizeof(gpc_vertex));
    gpc_vertex *p;
    
    for(i=n+1,p=v; --i; ++p) {
        p->x = *x++;
        p->y = *y++;
    }
    return v;
}

void mexFunction(int nL, mxArray *L[], int nR, const mxArray *R[])
{
    register int i;
    
    const int n1 = mxGetNumberOfElements(R[0]);
    const double *x1=mxGetPr(R[0]), *y1=mxGetPr(R[1]);
    
    const int n2 = mxGetNumberOfElements(R[2]);
    const double *x2=mxGetPr(R[2]), *y2=mxGetPr(R[3]);
    
    gpc_vertex_list l1 = {n1, mk_gpc_vertex(n1,x1,y1)};
    gpc_vertex_list l2 = {n2, mk_gpc_vertex(n2,x2,y2)};
    
    gpc_polygon p1 = {1, NULL, &l1}, p2 = {1, NULL, &l2}, p3;
    
    gpc_polygon_clip(GPC_INT, &p1, &p2, &p3);
    
    mwSize dims[] = {2, p3.contour[0].num_vertices};
    L[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    
    memcpy(mxGetPr(L[0]),p3.contour[0].vertex,2*dims[1]*sizeof(double));
    
    gpc_free_polygon(&p3);
    free(l2.vertex);
    free(l1.vertex);
}
