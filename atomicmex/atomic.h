#ifndef __ATOMIC_H__
#define __ATOMIC_H__

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct {
    const void *AP;
    void *TICS;
    void *diff_osc;
    double mu_2q;
} atomic_desc;

atomic_desc get_atomic_desc(const char *target, const char *type, const char *model);

double sigmav(double w, const atomic_desc *D);

void sigmav_vec(int n, double *w, const atomic_desc *D);


#ifdef __cplusplus
}
#endif

#endif /* __ATOMIC_H__ */
