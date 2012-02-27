#include <stdlib.h>

#include "tcp.c"

#include "mdsclient.h"

#include <ipdesc.h>

typedef struct {
    w_dtype_t w_dtype;
    char dtype;
} dtype_item;

static const dtype_item dtype_table[] = {
    {w_dtype_CSTRING       , DTYPE_CSTRING       },
    {w_dtype_UCHAR         , DTYPE_UCHAR         },
    {w_dtype_CHAR          , DTYPE_CHAR          },
    {w_dtype_USHORT        , DTYPE_USHORT        },
    {w_dtype_SHORT         , DTYPE_SHORT         },
    {w_dtype_ULONG         , DTYPE_ULONG         },
    {w_dtype_LONG          , DTYPE_LONG          },
    {w_dtype_ULONGLONG     , DTYPE_ULONGLONG     },
    {w_dtype_LONGLONG      , DTYPE_LONGLONG      },
    {w_dtype_FLOAT         , DTYPE_FLOAT         },
    {w_dtype_COMPLEX       , DTYPE_COMPLEX       },
    {w_dtype_DOUBLE        , DTYPE_DOUBLE        },
    {w_dtype_COMPLEX_DOUBLE, DTYPE_COMPLEX_DOUBLE}
};

static int dtype_table_len = sizeof(dtype_table)/sizeof(dtype_table[0]);


const char *get_dtype(w_dtype_t w_dtype)
{
    int i;
    for(i=0; i<dtype_table_len; i++)
        if (dtype_table[i].w_dtype == w_dtype) return &dtype_table[i].dtype;
    return NULL;
}

const w_dtype_t *get_w_dtype(char dtype)
{
    int i;
    for(i=0; i<dtype_table_len; i++)
        if (dtype_table[i].dtype == dtype) return &dtype_table[i].w_dtype;
    return NULL;
}


Descrip *mkDescrip(Descrip *l, w_dtype_t w_dtype, char ndims, int *dims, int num, int siz, void *ptr)
{
    int i;
    l->w_dtype = w_dtype;
    l->ndims = ndims;
    l->dims = (ndims==0) ? NULL : malloc(ndims*sizeof(int));
    for(i=0; i<ndims; i++) l->dims[i] = dims[i];
    l->num = num;
    l->siz = siz;
    l->ptr = ptr;
    return l;
}


int sm_mdsconnect(char *host) 
{
    return tcpconnect(host);
}

int sm_mdsdisconnect(int sock) 
{
    tcpdisconnect(sock);
    return 1;
}

int sm_mdsvalue(int sock, Descrip *l, int nr, Descrip *r, void **mem) 
{
    struct descrip exparg, *arg;
    const char *dtype;
    const w_dtype_t *w_dtype;
    int i, numbytes, num, siz, stat;

    for(i=0; i<nr; i++,r++) {
	dtype = get_dtype(r->w_dtype);
	arg = MakeDescrip(&exparg, *dtype, r->ndims, r->dims, r->ptr);
        stat = SendArg(sock, i, arg->dtype, nr, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    }

    *mem = NULL;
    stat = GetAnswerInfoTS(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr, mem);
    
    w_dtype = get_w_dtype(arg->dtype);
    siz = ArgLen(arg);
    num = numbytes/siz;

    mkDescrip(l, *w_dtype, arg->ndims, arg->dims, num, siz, arg->ptr);

    return 1;
}

int sm_mdsopen(int sock, char *tree, int shot) 
{
    return MdsOpen(sock,tree,shot);
}

int sm_mdsclose(int sock) 
{
    return MdsClose(sock);
}



