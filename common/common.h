#ifndef __COMMON__
#define __COMMON__

#define KV_LEN(x) (sizeof(x)/sizeof((x)[0]))

typedef double real;

typedef real (func)(real*);

typedef struct {
    const char *key;
    const void *val;
} keyval;

const void *kv_select(int n, const keyval *KV, const char *key);

#endif /* __COMMON__ */
