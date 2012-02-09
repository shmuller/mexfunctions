#include "stdlib.h"
#include "string.h"

#include "common.h"

const void *kv_select(int n, const keyval *KV, const char *key)
{
    int i;
    const keyval *kv;
    for(i=n,kv=KV; i--; kv++) {
        if (strcmp(key,kv->key)==0) {
            return kv->val;
        }
    }
    return NULL;
};
