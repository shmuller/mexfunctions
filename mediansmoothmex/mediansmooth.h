#ifdef __cplusplus
extern "C"
{
#endif

int median_filt_double(const double *x, const int N, const int w, int *ind);
int median_filt_float(const float *x, const int N, const int w, int *ind);

#ifdef __cplusplus
}
#endif
