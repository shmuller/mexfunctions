#ifndef __COMMON__
#define __COMMON__

#define LENGTH(x) (sizeof(x)/sizeof((x)[0]))

typedef double real;

typedef real (func)(real*);

typedef struct {
    const char *key;
    const void *val;
} keyval;

const void *select(int n, const keyval *KV, const char *key);

#endif /* __COMMON__ */
