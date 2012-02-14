#ifndef __COMMON_H__
#define __COMMON_H__

#ifdef __cplusplus
extern "C"
{
#endif

#define KV_LEN(x) (sizeof(x)/sizeof((x)[0]))

typedef double SM_REAL;

typedef SM_REAL (func)(SM_REAL*);

typedef struct {
    const char *key;
    const void *val;
} keyval;

const void *kv_select(int n, const keyval *KV, const char *key);

#ifdef __cplusplus
}  /* end extern "C" */
#endif

#endif /* __COMMON_H__ */
