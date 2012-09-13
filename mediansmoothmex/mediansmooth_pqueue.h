#ifdef __cplusplus
extern "C"
{
#endif

typedef int (cmpfun_t)(void *, void *);
typedef void (printfun_t)(void *);

typedef struct {
    cmpfun_t *lt;
    cmpfun_t *gt;
    printfun_t *print;
} fun_t;

void median_filt_pqueue(void *X, int N, int w, int *ind, int bdry, int bytes, fun_t *fun);

#ifdef __cplusplus
}
#endif
