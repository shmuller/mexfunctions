#ifndef IPDESC_H
#define IPDESC_H

#ifdef _WIN32
#ifndef _WS2DEF_
#include "windows.h"
#endif
#else
typedef int SOCKET;
#endif

#define MAX_DIMS 7
#define DTYPE_UCHAR   2
#define DTYPE_USHORT  3
#define DTYPE_ULONG   4
#define DTYPE_ULONGLONG 5
#define DTYPE_CHAR    6
#define DTYPE_SHORT   7
#define DTYPE_LONG    8
#define DTYPE_LONGLONG 9
#ifdef DTYPE_FLOAT
#undef DTYPE_FLOAT
#endif
#define DTYPE_FLOAT   10
#ifdef DTYPE_DOUBLE
#undef DTYPE_DOUBLE
#endif
#define DTYPE_DOUBLE  11
#define DTYPE_COMPLEX 12
#define DTYPE_COMPLEX_DOUBLE 13
#define DTYPE_CSTRING 14
#define DTYPE_EVENT_NOTIFY   99
#ifndef DTYPE_EVENT
#define DTYPE_EVENT DTYPE_EVENT_NOTIFY
#endif

struct descrip {
    char dtype;
    char ndims;
    int  dims[MAX_DIMS];
    short length;
	void *ptr;
};

#ifdef __cplusplus
extern "C" {
#endif

struct descrip *MakeDescrip(struct descrip *in_descrip, char dtype, 
    char ndims, int *dims, void *ptr);

short ArgLen(struct descrip *d);

int SendArg(SOCKET sock, unsigned char idx, char dtype, unsigned char nargs, 
    short length, char ndims, int *dims, char *bytes);

int GetAnswerData(SOCKET sock, char *dtype, short *length, char *ndims, int *dims, 
    int *numbytes, void **dptr);


#ifdef __cplusplus
} // extern "C"
#endif

#endif
