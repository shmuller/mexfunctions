#ifndef __MDSCLIENT_H__
#define __MDSCLIENT_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    w_dtype_UNKNOWN=-1,
    w_dtype_CSTRING,
    w_dtype_UCHAR,
    w_dtype_CHAR,
    w_dtype_USHORT,
    w_dtype_SHORT,
    w_dtype_ULONG,
    w_dtype_LONG,
    w_dtype_ULONGLONG,
    w_dtype_LONGLONG,
    w_dtype_FLOAT,
    w_dtype_COMPLEX,
    w_dtype_DOUBLE,
    w_dtype_COMPLEX_DOUBLE
} w_dtype_t;

typedef struct {
    w_dtype_t w_dtype;
    char ndims;
    int *dims;
    int num;
    int siz;
    void *ptr;
} Descrip;

Descrip *mkDescrip_dims(Descrip *l, char ndims, int *dims, int num, int siz);
Descrip *mkDescrip_data(Descrip *l, w_dtype_t w_dtype, void *ptr);
Descrip *mkDescrip(Descrip *l, w_dtype_t w_dtype, char ndims, int *dims, int num, int siz, void *ptr);

int sm_mdsvalue(int sock, Descrip *l, int nr, Descrip *r, void **mem);

int sm_mdsconnect(char *host);
int sm_mdsdisconnect(int sock);

#ifdef __cplusplus
}
#endif

#endif /* __MDSCLIENT_H__ */

