#ifndef ___TRIBLOCK_H__
#define ___TRIBLOCK_H__

#ifdef __cplusplus
extern "C" {
#endif

void _triblock(double *D1, double *U2, int *IPIV, double *L1, double *D2, int is_half,
    int N, int NRHS, int *INFO);

#ifdef __cplusplus
}
#endif

#endif /* ___TRIBLOCK_H__ */

