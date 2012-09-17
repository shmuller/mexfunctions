#ifdef __cplusplus
extern "C"
{
#endif

typedef int (cmpfun_t)(void *, void *);
typedef void (printfun_t)(void *);
typedef double (to_double_t)(void *);
typedef void (from_double_t)(void *, double);

typedef struct {
    cmpfun_t *lt;
    cmpfun_t *gt;
    printfun_t *print;
    to_double_t *to_double;
    from_double_t *from_double;
} fun_t;

void median_filt_pqueue(void *X, int N, int w, int *ind, int bdry, int bytes, fun_t *fun);

#ifdef __cplusplus
}
#endif
