#include <stdlib.h>
#include <string.h>

#include "tcp.h"

#include "mdsclient.h"

#include "ipdesc.h"

#ifndef status_ok
#define status_ok(status) (((status) & 1) == 1)
#endif


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


char get_dtype(w_dtype_t w_dtype)
{
    int i;
    for(i=0; i<dtype_table_len; i++)
        if (dtype_table[i].w_dtype == w_dtype) return dtype_table[i].dtype;
    return 0;
}

w_dtype_t get_w_dtype(char dtype)
{
    int i;
    for(i=0; i<dtype_table_len; i++)
        if (dtype_table[i].dtype == dtype) return dtype_table[i].w_dtype;
    return w_dtype_UNKNOWN;
}


Descrip *mkDescrip_dims(Descrip *l, char ndims, int *dims, int num, int siz)
{
    int i;
    l->ndims = ndims;
    l->dims = (ndims==0) ? NULL : malloc(ndims*sizeof(int));
    for(i=0; i<ndims; i++) l->dims[i] = dims[i];
    l->num = num;
    l->siz = siz;
    return l;
}

Descrip *mkDescrip_data(Descrip *l, w_dtype_t w_dtype, void *ptr)
{
    l->w_dtype = w_dtype;
    l->ptr = ptr;
    return l;
}

Descrip *mkDescrip(Descrip *l, w_dtype_t w_dtype, char ndims, int *dims, int num, int siz, void *ptr)
{
    mkDescrip_dims(l, ndims, dims, num, siz);
    mkDescrip_data(l, w_dtype, ptr);
    return l;
}


int sm_mdsvalue(int sock, Descrip *l, int nr, Descrip *r, void **mem) 
{
    struct descrip exparg, *arg;
    char dtype;
    w_dtype_t w_dtype;
    int i, numbytes, num, siz, stat;

    for(i=0; i<nr; i++,r++) {
        dtype = get_dtype(r->w_dtype);
        arg = MakeDescrip(&exparg, dtype, r->ndims, r->dims, r->ptr);
        stat = SendArg(sock, i, arg->dtype, nr, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    }

    stat = GetAnswerData(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr);

    w_dtype = get_w_dtype(arg->dtype);
    siz = (w_dtype==w_dtype_CSTRING) ? sizeof(char) : ArgLen(arg);
    num = (siz==0) ? 0 : numbytes/siz;

    mkDescrip(l, w_dtype, arg->ndims, arg->dims, num, siz, arg->ptr);
    *mem = l->ptr;
    return 1;
}


int tcpauth(int sock, char *user_p)
{  
    struct descrip exparg, *arg;
    int numbytes = 0, stat = 0;
	
    arg = MakeDescrip(&exparg,DTYPE_CSTRING,0,NULL,user_p);
    stat = SendArg(sock, 0, arg->dtype, 1, ArgLen(arg), arg->ndims, arg->dims, arg->ptr);
    stat = GetAnswerData(sock, &arg->dtype, &arg->length, &arg->ndims, arg->dims, &numbytes, &arg->ptr);
    if (arg->ptr) free(arg->ptr);
    if (!status_ok(stat)) {
        tcpclose(sock);
        return -5;
    }
    return 1;
}


int sm_mdsconnect(char *host) 
{
    char *port;
    int sock, err;
    int STRLEN = 4096;
    char user[STRLEN];
    char *user_p = host;

    if ((host=strchr(user_p,'@')) == NULL) {
        host = user_p;
        user_p = tcpuser(user,STRLEN);
    } else {
        *host++ = 0;
    }
    if ((port=strchr(host,':')) == NULL) {
        port = strdup("8000");
    } else {
        *port++ = 0;
    }

    if ((sock=tcpopen(host,port)) < 0) {
        return sock;		
    }
    if ((err=tcpauth(sock,user_p)) < 0) {
        return err;
    }
    return sock;

}

int sm_mdsdisconnect(int sock) 
{
    tcpclose(sock);
    return 1;
}


