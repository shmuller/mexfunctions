#ifdef __cplusplus
extern "C"
{
#endif

typedef int (cmpfun_t)(void *, void *);
typedef void (printfun_t)(void *);
typedef double (as_double_t)(void *);

typedef struct {
    cmpfun_t *lt;
    cmpfun_t *gt;
    printfun_t *print;
    as_double_t *as_double;
} fun_t;

void median_filt_pqueue(void *X, int N, int w, int *ind, int bdry, int bytes, fun_t *fun);

#ifdef __cplusplus
}
#endif
